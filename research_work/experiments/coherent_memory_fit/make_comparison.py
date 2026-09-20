"""Export the requested observed/candidate/Newtonian comparisons."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import gradient_laws as L
D=L.D;HERE=L.HERE
def main():
    sel=D.read(HERE/"evidence-gradient-v1"/"selection-before-test.json")["selected"]["joint"]
    force=L.candidate_force(sel)
    # Named example selected before inspecting its fit: a well measured extended disk.
    g=next(g for g in D.galaxies() if g["name"]=="NGC3198")
    vnew=np.sqrt(g["r"]*force(g));vnewton=np.sqrt(g["r"]*g["gb"])
    observed=D.coma_data()
    radii=np.linspace(min(observed["r"])*.85,max(observed["r"])*1.03,180)
    models=[];source_metadata=[]
    for ne0,mstar in ((2.5e-3,.5e13),(4.5e-3,2e13)):
        src=D.coma_source(ne0,mstar,8192)
        law=L.optical_law(src,force,9000,8192)
        standard=L.optical_law(src,lambda o:o["gb"],9000,8192)
        cand=L.OLD.shear_shape(law,radii,128,.002,(3000.,9000.))["shape"]*50000
        ordinary=L.OLD.shear_shape(standard,radii,128,.002,(3000.,9000.))["shape"]*50000
        models.append(dict(candidate=cand,newtonian_ballistic=.5*ordinary,ordinary_standard_light=ordinary))
        source_metadata.append(dict(ne0_cm3=ne0,stellar_mass_Msun=mstar,total_ordinary_mass_Msun=src["mass"]))
    fig,axes=plt.subplots(1,2,figsize=(15,7.4),gridspec_kw={"wspace":.25})
    fig.patch.set_facecolor("#f7f9fc")
    blue="#1f66c1";orange="#c56725";gray="#8492a2"
    for ax in axes:
        ax.set_facecolor("white");ax.spines[["top","right"]].set_visible(False)
        ax.grid(alpha=.18);ax.tick_params(labelsize=10)
    ax=axes[0];x=radii/1000
    for key,color,label,alpha in [
        ("candidate",blue,"Our shared candidate: matter bracket",.2),
        ("newtonian_ballistic",orange,"Newtonian ballistic light: matter bracket",.2)]:
        ys=np.array([m[key] for m in models])*1000
        ax.fill_between(x,np.min(ys,axis=0),np.max(ys,axis=0),color=color,alpha=alpha,label=label)
        for y in ys:ax.plot(x,y,color=color,linewidth=1.5,alpha=.9)
    for i,m in enumerate(models):
        ax.plot(x,m["ordinary_standard_light"]*1000,":",color=gray,linewidth=1.8,
            label="Ordinary matter + standard light bending" if i==0 else None)
    ax.errorbar(observed["r"]/1000,observed["y"]*1000,yerr=observed["error"]*1000,
        fmt="o",color="#20242c",capsize=4,markersize=5,label="Observed shear ± plotted error",zorder=8)
    ax.axhline(0,color="#444",linewidth=.8)
    ax.set_title("Coma cluster | light-image distortion",loc="left",fontsize=14,fontweight="bold",pad=34)
    ax.text(0,1.02,"Fixed 100 Mpc lens distance; source efficiency beta = 1",transform=ax.transAxes,fontsize=10,color="#516173")
    ax.set_xlabel("Projected distance from cluster center (Mpc)",labelpad=9)
    ax.set_ylabel("Tangential shear (× 0.001)",labelpad=9)
    ax.legend(loc="upper right",frameon=True,fontsize=8.2)
    ax.set_xlim(x[0],x[-1]);ax.xaxis.set_major_locator(MaxNLocator(7))
    ax=axes[1]
    ax.plot(g["r"],vnew,color=blue,lw=2.6,label="Our same shared candidate",zorder=4)
    ax.plot(g["r"],vnewton,color=orange,lw=2.6,label="Newtonian, ordinary matter only")
    ax.errorbar(g["r"],g["y"],yerr=g["error"],fmt="o",ms=3.4,color="#20242c",
        capsize=2,label="Observed rotation ± quoted error",zorder=5)
    ax.axvspan(g["r"][-1]*.75,g["r"][-1],color=blue,alpha=.06,label="Outer measured disk")
    ax.set_title("NGC 3198 | circular orbital speed",loc="left",fontsize=14,fontweight="bold",pad=34)
    ax.text(0,1.02,f"Fixed catalog distance {g['catalog']['distance']:g} Mpc; {g['split']}-labelled object",
        transform=ax.transAxes,fontsize=10,color="#516173")
    ax.set_xlabel("Distance from galaxy center (kpc)",labelpad=9)
    ax.set_ylabel("Circular speed (km/s)",labelpad=9)
    ax.set_ylim(bottom=0);ax.set_xlim(0,g["r"][-1]*1.04)
    ax.legend(loc="lower right",fontsize=8.5)
    fig.suptitle("One shared formula, compared with measurements",x=.06,ha="left",y=.97,fontsize=20,fontweight="bold")
    fig.text(.06,.895,"Candidate CMF-1/2: RM-ell0-w1-ridge0. No object-specific refit, dark matter, expansion or adjusted distances.",
        fontsize=10.5,color="#445469")
    fig.subplots_adjust(left=.06,right=.98,bottom=.28,top=.78)
    fig.text(.06,.13,"Coma: bands are two assumed ordinary-source profiles, not confidence intervals. beta = 1 is the maximum static source-distance factor;\n"
        "finite source distances reduce all shown lensing curves. No fitted lens amplitude is used. Points were reconstructed from Kubo et al.'s figure.",
        fontsize=9.2,color="#445469",linespacing=1.5)
    fig.text(.06,.06,"Newtonian light baseline: ballistic rays at speed c, giving half the standard weak-field deflection. The dotted gray curves show that standard\n"
        "ordinary-matter light law for context. Galaxy points trace the rotation field (often gas), not an individually tracked star. All data were previously exposed.",
        fontsize=9.2,color="#445469",linespacing=1.5)
    for suffix in ("png","svg","pdf"):
        fig.savefig(HERE/("observed-comparison."+suffix),dpi=170,facecolor=fig.get_facecolor())
    svg=HERE/"observed-comparison.svg"
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text(encoding="utf8").splitlines())+"\n",encoding="utf8",newline="\n")
    plt.close(fig)
    L.save(HERE/"observed-comparison-data.json",dict(candidate=sel,galaxy_selection="NGC3198 named before inspecting fit; illustrative, not a statistical winner",
        galaxy=dict(name=g["name"],split=g["split"],distance_Mpc=g["catalog"]["distance"],radius_kpc=g["r"],
            observed_kms=g["y"],error_kms=g["error"],candidate_kms=vnew,newtonian_kms=vnewton),
        coma=dict(data=observed,radius_kpc=radii,models=models,source_brackets=source_metadata,distance_Mpc=100,beta=1,
            geometry="Conditional static Euclidean geometry; archived radius conversion remains approximate.",
            band_interpretation="ordinary-source bracket, not confidence interval",
            newtonian_light="ballistic v=c test ray; half the stipulated standard weak-field light deflection")))
    print(json.dumps(dict(galaxy=g["name"],split=g["split"],last_radius_kpc=float(g["r"][-1]),
        observed=float(g["y"][-1]),error=float(g["error"][-1]),candidate=float(vnew[-1]),newtonian=float(vnewton[-1])),indent=2))
if __name__=="__main__":main()
