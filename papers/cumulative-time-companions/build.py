"""Build the working manuscript PDF from Markdown; no research outputs are changed."""
from pathlib import Path
import argparse
import html
import io
import re

import matplotlib
matplotlib.use("Agg")
from matplotlib import mathtext
from matplotlib.font_manager import FontProperties
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.utils import ImageReader
from reportlab.platypus import SimpleDocTemplate, Paragraph, Flowable, PageBreak, KeepTogether, Image, Table, TableStyle


class Equation(Flowable):
    def __init__(self, expression):
        super().__init__()
        label = re.search(r"\\tag\{(\d+)\}", expression)
        self.label = label.group(1) if label else ""
        expression = re.sub(r"\\tag\{\d+\}", "", expression).strip()
        self.buffer = io.BytesIO()
        mathtext.math_to_image("$" + expression + "$", self.buffer, dpi=240, format="png",
                              prop=FontProperties(size=12))
        self.buffer.seek(0)
        self.picture = ImageReader(self.buffer)
        self.iw, self.ih = self.picture.getSize()

    def wrap(self, available_width, available_height):
        self.width = available_width
        self.draw_width = min(self.iw * 72 / 240, available_width - 38)
        self.draw_height = self.draw_width * self.ih / self.iw
        self.height = self.draw_height + 23
        return self.width, self.height

    def draw(self):
        self.canv.drawImage(self.picture, (self.width - self.draw_width) / 2 - 10,
                            11, self.draw_width, self.draw_height, mask="auto")
        self.canv.setFont("Times-Roman", 10)
        self.canv.drawRightString(self.width, self.height / 2, f"({self.label})")


def rich(text):
    text = html.escape(text)
    text = re.sub(r"\[([^\]]+)\]\((https?://[^)]+)\)",
                  r'<link href="\2" color="#234f77">\1</link>', text)
    return text


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("manuscript.pdf"))
    args = parser.parse_args()
    text = Path(__file__).with_name("manuscript.md").read_text(encoding="utf-8")
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle("PaperBody", fontName="Times-Roman", fontSize=11,
                              leading=15, spaceAfter=8))
    styles.add(ParagraphStyle("PaperTitle", fontName="Times-Bold", fontSize=24,
                              leading=28, spaceAfter=14, textColor=colors.HexColor("#16334a")))
    styles.add(ParagraphStyle("PaperSubtitle", fontName="Times-Roman", fontSize=15,
                              leading=19, spaceAfter=14))
    styles.add(ParagraphStyle("PaperSection", fontName="Times-Bold", fontSize=13,
                              leading=17, spaceBefore=14, spaceAfter=7, keepWithNext=True))
    styles.add(ParagraphStyle("PaperSubsection", fontName="Times-Bold", fontSize=11,
                              leading=15, spaceBefore=9, spaceAfter=6, keepWithNext=True))
    styles.add(ParagraphStyle("PaperBullet", parent=styles["PaperBody"], leftIndent=13,
                              firstLineIndent=-9, spaceAfter=5))
    story = []
    equation_count = 0
    blocks = re.split(r"\n\s*\n", text.strip())
    for index, block in enumerate(blocks):
        if block.startswith("$$"):
            equation_count += 1
            story.append(Equation(block.strip().removeprefix("$$").removesuffix("$$").strip()))
        elif block.startswith("!["):
            match = re.fullmatch(r'!\[([^\]]*)\]\(([^)]+)\)', block)
            assert match, block
            picture = Image(str(Path(__file__).parent / match.group(2)))
            picture.drawHeight *= 496 / picture.drawWidth
            picture.drawWidth = 496
            story.append(KeepTogether([picture, Paragraph(rich(match.group(1)), styles['PaperBody'])]))
        elif block.startswith('|'):
            rows = [[Paragraph(rich(cell.strip()), styles['PaperBullet']) for cell in line.strip('|').split('|')]
                    for line in block.splitlines() if not re.match(r'^\|[\s:|\-]+\|$',line)]
            table = Table(rows, repeatRows=1, hAlign='LEFT', colWidths=[496/len(rows[0])]*len(rows[0]))
            table.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e4edf3')),
                ('VALIGN',(0,0),(-1,-1),'TOP'),('LINEBELOW',(0,0),(-1,0),.7,colors.grey),
                ('BOTTOMPADDING',(0,0),(-1,-1),7),('TOPPADDING',(0,0),(-1,-1),7)]))
            table.keepWithNext = True
            story.append(KeepTogether([table]))
        elif block.startswith("#### "):
            story.append(Paragraph(rich(block[5:]), styles["PaperSubsection"]))
        elif block.startswith("### "):
            if block == "### References":
                story.append(PageBreak())
            story.append(Paragraph(rich(block[4:]), styles["PaperSection"]))
        elif block.startswith("## "):
            story.append(Paragraph(rich(block[3:]), styles["PaperSubtitle"]))
        elif block.startswith("# "):
            story.append(Paragraph(rich(block[2:]), styles["PaperTitle"]))
        elif block.startswith("- "):
            for line in block.splitlines():
                story.append(Paragraph("&#8226; " + rich(line[2:]), styles["PaperBullet"]))
        else:
            paragraph = Paragraph(rich(" ".join(block.splitlines())), styles["PaperBody"])
            if index + 1 < len(blocks) and blocks[index + 1].startswith("$$"):
                paragraph.keepWithNext = True
            story.append(KeepTogether([paragraph]) if re.match(r"\[\d\] ", block) else paragraph)
    assert equation_count == 20
    args.output.parent.mkdir(parents=True, exist_ok=True)
    doc = SimpleDocTemplate(str(args.output), pagesize=(612, 792), rightMargin=58,
                            leftMargin=58, topMargin=49, bottomMargin=48,
                            title="Environmental Time Stretching and Companion-Energy Deposition",
                            author="", subject="Working theoretical framework, version 0.5")

    def page(canvas, document):
        canvas.saveState()
        canvas.setFont("Times-Roman", 8)
        canvas.setFillColor(colors.HexColor("#526371"))
        canvas.drawString(58, 767, "ENVIRONMENTAL TIME STRETCHING AND COMPANION ENERGY")
        canvas.drawString(58, 28, "Working draft v0.5  |  10 September 2026  |  Theory incomplete")
        canvas.drawRightString(554, 28, str(document.page))
        canvas.restoreState()

    doc.build(story, onFirstPage=page, onLaterPages=page)
    print(args.output)


if __name__ == "__main__":
    main()
