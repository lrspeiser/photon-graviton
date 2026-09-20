"""Render the project research brief from its editable JSON source."""
from pathlib import Path
import json,html,re,hashlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,PageBreak,Table,TableStyle,Image,Preformatted,KeepTogether
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from pypdf import PdfReader
ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/"research_plan/companion-stream-research-brief.json"
OUT=ROOT/"output/pdf/companion-stream-gravity-research-brief.pdf"
TMP=ROOT/"tmp/pdfs/companion-stream-brief"
EXP=ROOT/"research_work/experiments/coherent_memory_fit"
NAVY=colors.HexColor("#153349");TEAL=colors.HexColor("#176c80")
def esc(s):return html.escape(str(s)).replace("\n","<br/>")
def figures():
    d=json.loads((EXP/"observed-comparison-data.json").read_text(encoding="utf8"))
    fig,axs=plt.subplots(1,2,figsize=(7.2,3.35))
    blue="#246ac0";orange="#c56a2c";gray="#8292a1"
    c=d["coma"];x=np.array(c["radius_kpc"])/1000
    for key,col in (("candidate",blue),("newtonian_ballistic",orange)):
        arr=np.array([m[key] for m in c["models"]])*1000
        axs[0].fill_between(x,arr.min(axis=0),arr.max(axis=0),color=col,alpha=.23)
        for y in arr:axs[0].plot(x,y,color=col,lw=1.2)
    for m in c["models"]:
        axs[0].plot(x,np.array(m["ordinary_standard_light"])*1000,":",color=gray,lw=1.2)
    o=c["data"];axs[0].errorbar(np.array(o["r"])/1000,np.array(o["y"])*1000,
        yerr=np.array(o["error"])*1000,fmt="o",ms=3,color="#20262d",capsize=2)
    axs[0].axhline(0,c="#888",lw=.6);axs[0].set_title("Coma: lensing shear",loc="left",fontsize=9.5,fontweight="bold")
    axs[0].set_xlabel("Projected radius (Mpc)",fontsize=8);axs[0].set_ylabel("Shear (x 0.001)",fontsize=8)
    g=d["galaxy"]
    axs[1].plot(g["radius_kpc"],g["candidate_kms"],color=blue,lw=1.8)
    axs[1].plot(g["radius_kpc"],g["newtonian_kms"],color=orange,lw=1.8)
    axs[1].errorbar(g["radius_kpc"],g["observed_kms"],yerr=g["error_kms"],fmt="o",color="#20262d",ms=2,capsize=1)
    axs[1].set_title("NGC 3198: rotation",loc="left",fontsize=9.5,fontweight="bold")
    axs[1].set_xlabel("Galactic radius (kpc)",fontsize=8);axs[1].set_ylabel("Circular speed (km/s)",fontsize=8)
    axs[1].set_ylim(bottom=0)
    for ax in axs:
        ax.spines[["right","top"]].set_visible(False);ax.grid(alpha=.16);ax.tick_params(labelsize=7.5)
    fig.text(.04,.945,"Blue: shared candidate   |   Black: observed   |   Orange: Newtonian",fontsize=8.6,color="#31485a")
    fig.text(.04,.015,"Coma dotted gray: standard light bending from ordinary matter. No dark matter in any plotted curve.",
        fontsize=7.3,color="#516173")
    fig.subplots_adjust(left=.09,right=.99,top=.82,bottom=.18,wspace=.4)
    fig.savefig(TMP/"comparison.png",dpi=240)
    plt.close(fig)
def markdown(data):
    lines=["# "+data["title"],"",data["subtitle"],"",data["date"]+" | Fictional-universe research programme","",
           "Evidence snapshot: "+data["evidence_commit"],""]
    for page in data["pages"]:
        lines+=["## "+page["title"],""]
        for b in page["blocks"]:
            kind=b["type"]
            if kind=="h":lines+=["### "+b["text"],""]
            elif kind=="equation":lines+=["\x60\x60\x60text",b["text"],"\x60\x60\x60",""]
            elif kind=="table":
                lines+=["| "+" | ".join(s.replace("\n"," ") for s in b["headers"])+" |",
                        "| "+" | ".join("---" for _ in b["headers"])+" |"]
                lines+=["| "+" | ".join(r)+" |" for r in b["rows"]];lines+=[""]
            elif kind=="figure":lines+=["![Observed comparison](../research_work/experiments/coherent_memory_fit/observed-comparison.png)",""]
            elif kind=="source":lines+=[f"[{b['id']}] {b['text']} [Read source]({b['url']}).",""]
            elif kind=="callout":lines+=["> "+b["text"],""]
            else:lines+=[b["text"],""]
    (ROOT/"research_plan/companion-stream-research-brief.md").write_text("\n".join(lines),encoding="utf8",newline="\n")
def main():
    data=json.loads(SOURCE.read_text(encoding="utf8"));TMP.mkdir(parents=True,exist_ok=True);OUT.parent.mkdir(parents=True,exist_ok=True)
    figures();markdown(data)
    styles=getSampleStyleSheet()
    styles.add(ParagraphStyle(name="BodyBrief",fontName="Helvetica",fontSize=10.1,leading=14.1,textColor=NAVY,spaceAfter=8))
    styles.add(ParagraphStyle(name="LeadBrief",parent=styles["BodyBrief"],fontSize=11.5,leading=15.7,textColor=TEAL,spaceAfter=12))
    styles.add(ParagraphStyle(name="HBrief",parent=styles["BodyBrief"],fontName="Helvetica-Bold",fontSize=11.3,leading=14,spaceBefore=6,spaceAfter=6,keepWithNext=True))
    styles.add(ParagraphStyle(name="TitleBrief",fontName="Helvetica-Bold",fontSize=23,leading=27,textColor=NAVY,spaceAfter=12,keepWithNext=True))
    styles.add(ParagraphStyle(name="SmallBrief",parent=styles["BodyBrief"],fontSize=8.1,leading=10.7,textColor=colors.HexColor("#4d6274"),spaceAfter=7))
    styles.add(ParagraphStyle(name="SourceBrief",parent=styles["SmallBrief"],fontSize=8.8,leading=11.4,spaceAfter=9))
    styles.add(ParagraphStyle(name="CellBrief",parent=styles["SmallBrief"],fontSize=8.7,leading=11,textColor=NAVY,spaceAfter=0))
    styles.add(ParagraphStyle(name="HeadCell",parent=styles["CellBrief"],fontName="Helvetica-Bold",textColor=colors.white))
    styles.add(ParagraphStyle(name="EquationBrief",fontName="Courier",fontSize=9,leading=12,textColor=NAVY,spaceBefore=3,spaceAfter=10))
    story=[]
    for i,page in enumerate(data["pages"]):
        if i:story.append(PageBreak())
        story.append(Paragraph(f"RESEARCH BRIEF / {i+1:02d}",styles["SmallBrief"]))
        story.append(Paragraph(esc(page["title"]),styles["TitleBrief"]))
        if i==0:
            story.append(Paragraph(esc(data["title"]+" - "+data["subtitle"]),styles["HBrief"]))
            story.append(Paragraph("Fictional-universe research programme | "+data["date"],styles["SmallBrief"]))
        for b in page["blocks"]:
            kind=b["type"]
            if kind in ("p","lead","h","small"):
                sty={"p":"BodyBrief","lead":"LeadBrief","h":"HBrief","small":"SmallBrief"}[kind]
                story.append(Paragraph(esc(b["text"]),styles[sty]))
            elif kind=="equation":
                story.append(Preformatted(b["text"],styles["EquationBrief"]))
            elif kind=="callout":
                t=Table([[Paragraph(esc(b["text"]),styles["BodyBrief"])]],colWidths=[499])
                t.setStyle(TableStyle([("BACKGROUND",(0,0),(-1,-1),colors.HexColor("#edf5f8")),
                    ("BOX",(0,0),(-1,-1),.5,colors.HexColor("#bad7e1")),
                    ("LEFTPADDING",(0,0),(-1,-1),10),("RIGHTPADDING",(0,0),(-1,-1),10),
                    ("TOPPADDING",(0,0),(-1,-1),9),("BOTTOMPADDING",(0,0),(-1,-1),3)]))
                story.extend([t,Spacer(1,10)])
            elif kind=="table":
                table=[[Paragraph(esc(v),styles["HeadCell"]) for v in b["headers"]]]
                table += [[Paragraph(esc(v),styles["CellBrief"]) for v in row] for row in b["rows"]]
                t=Table(table,colWidths=b["widths"],repeatRows=1,hAlign="LEFT")
                t.setStyle(TableStyle([("BACKGROUND",(0,0),(-1,0),TEAL),
                    ("ROWBACKGROUNDS",(0,1),(-1,-1),[colors.HexColor("#f0f5f8"),colors.white]),
                    ("VALIGN",(0,0),(-1,-1),"TOP"),("LEFTPADDING",(0,0),(-1,-1),8),
                    ("TOPPADDING",(0,0),(-1,-1),6),("BOTTOMPADDING",(0,0),(-1,-1),6),
                    ("LINEBELOW",(0,-1),(-1,-1),.5,colors.HexColor("#b9c9d5"))]))
                story.extend([t,Spacer(1,10)])
            elif kind=="figure":
                story.append(Image(str(ROOT/b["path"]),width=499,height=499*3.35/7.2))
                story.append(Spacer(1,6))
            elif kind=="source":
                text=f'<b>[{b["id"]}]</b> '+esc(b["text"])+f' <link href="{html.escape(b["url"],quote=True)}" color="#176c80">Read source</link>.'
                story.append(Paragraph(text,styles["SourceBrief"]))
    def footer(canvas,doc):
        canvas.saveState();w,h=A4
        canvas.setStrokeColor(colors.HexColor("#c6d3db"));canvas.line(48,43,w-48,43)
        canvas.setFont("Helvetica",8);canvas.setFillColor(colors.HexColor("#516173"))
        canvas.drawString(48,30,"PHOTON-GRAVITON  |  Research proposal  |  Evidence: f62f518")
        canvas.drawRightString(w-48,30,str(doc.page))
        canvas.restoreState()
    doc=SimpleDocTemplate(str(OUT),pagesize=A4,rightMargin=48,leftMargin=48,topMargin=42,bottomMargin=57,
        title=data["title"]+": "+data["subtitle"],author="Photon-graviton project",subject="Research status and proposed next tests")
    doc.build(story,onFirstPage=footer,onLaterPages=footer)
    reader=PdfReader(OUT);texts=[p.extract_text() or "" for p in reader.pages]
    if len(reader.pages)!=len(data["pages"]):raise AssertionError(f"Layout overflow: expected 7 pages, got {len(reader.pages)}")
    for text,page in zip(texts,data["pages"]):
        if page["title"] not in text:raise AssertionError("Missing or displaced page title: "+page["title"])
    for marker in ("149","119.3","64.7","5.06","58.06","2.6 sigma","Priority 6"):
        if marker not in "\n".join(texts):raise AssertionError("Missing expected brief text: "+marker)
    digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    report=dict(pages=len(reader.pages),source_sha256=digest(SOURCE),pdf_sha256=digest(OUT),
        evidence_commit=data["evidence_commit"],evidence_files={str(p.relative_to(ROOT)):digest(p) for p in [
            EXP/"summary.json",EXP/"observed-comparison-data.json",EXP/"audit.json"]},
        page_text_lengths=[len(t) for t in texts],all_text_checks_passed=True,
        visual_review_required=True)
    (TMP/"text-and-layout-checks.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf8")
    print(json.dumps(report,indent=2))
if __name__=="__main__":main()
