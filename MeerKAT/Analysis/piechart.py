import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

total = [80,15,5]
labels = ["Dark Matter", "Intracluster Medium", "Galaxies"]
plt.title('Galaxy cluster mass budget')
plt.gca().axis("equal")
wedges, texts = plt.pie(total, startangle=90, labels=labels,
                        wedgeprops = { 'linewidth': 2, "edgecolor" :"k","fill":False,  })


def img_to_pie( fn, wedge, xy, zoom=1, ax = None):
    if ax==None: ax=plt.gca()
    im = plt.imread(fn, format='jpg')
    path = wedge.get_path()
    patch = PathPatch(path, facecolor='none')
    ax.add_patch(patch)
    imagebox = OffsetImage(im, zoom=zoom, clip_path=patch, zorder=-10)
    ab = AnnotationBbox(imagebox, xy, xycoords='data', pad=0, frameon=False)
    ax.add_artist(ab)

positions = [(-1,0.3),(0,-0.5),(0.5,0.5)]
zooms = [0.4,0.4,0.4]
files = ["univers_ce8bef718e.jpg", "univers_ce8bef718e.jpg","Bullet_cluster-2.jpg"]
for i in range(3):
    fn = "/net/rijn9/data2/swart/DATA/MeerKAT_DATA/Analysis/DATA/{}".format(files[i].lower())
    img_to_pie(fn, wedges[i], xy=positions[i], zoom=zooms[i] )
    wedges[i].set_zorder(10)

plt.show()
