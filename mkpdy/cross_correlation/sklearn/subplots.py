fig, axes = plt.subplots(nrows=2, ncols=2)
for ax in axes.flat:
    im = ax.imshow(emp_cov_residue, interpolation='gaussian', vmin=-1, vmax=1,cmap=plt.cm.RdBu_r, aspect="auto")
    #fig.colorbar(im)

axes.flat[0].set_xlim(520, 620)
axes.flat[0].set_ylim(700, 800)
axes.flat[1].set_xlim(450, 540)
axes.flat[1].set_ylim(900, 1000)
axes.flat[2].set_xlim(980, 1080)
axes.flat[2].set_ylim(710, 810)
axes.flat[3].set_xlim(520, 600)
axes.flat[3].set_ylim(850, 900)



plt.tight_layout()
#plt.savefig("covar.pdf")
plt.show()