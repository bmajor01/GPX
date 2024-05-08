import gpxpy
from gpxslicer import slicer
import matplotlib.pyplot as plt
import numpy as np
from gpxplotter import read_gpx_file, create_folium_map, add_segment_to_map
import folium
from math import pi, sqrt, exp
import matplotlib as mpl

from matplotlib.collections import LineCollection
from matplotlib.colors import BoundaryNorm, ListedColormap


def gauss(n,sigma):
    r = range(-int(n/2),int(n/2)+1)
    return [1 / (sigma * sqrt(2*pi)) * exp(-float(x)**2/(2*sigma**2)) for x in r]

#kernel = gauss(15,2)
n = 2
kernel = np.ones(n)
kernel = kernel / n

m = folium.Map(location=(47.1625, 19.5033),zoom_start=8)
folium.TileLayer('openstreetmap', name='OpenStreet Map').add_to(m)

gpx_file = open('okt_teljes_20230921.gpx')
gpx = gpxpy.parse(gpx_file)

print("{} track(s)".format(len(gpx.tracks)))
track = gpx.tracks[0]

print("{} segment(s)".format(len(track.segments)))
segment = track.segments[0]

print("{} point(s)".format(len(segment.points)))

print(track.length_3d())

days = 1

sliced = slicer.slice_gpx_at_interval(gpx, (track.length_3d()/days))

prevElevation = 0
prevFiltElevation = 0

fname = 1

for n in sliced.tracks:
    elevgain = 0
    filtElevgain = 0
    x = []
    segment = n.segments[0]

    points = []
    step = 1
    for point in segment.points[::step]:
        points.append(tuple([point.latitude, point.longitude]))

    folium_gpx = folium.PolyLine(points, color='red', weight=5, opacity=0.85).add_to(m)
    folium.TileLayer('openstreetmap', name='OpenStreet Map').add_to(m)
    m.save('KK_day'+str(fname)+'.html')
    import webbrowser

    #webbrowser.open('KK_day'+str(fname)+'.html')

    m = folium.Map(location=(47.1625, 19.5033),zoom_start=8)

    startPoint = segment.points[0]
    prevElevation = startPoint.elevation

    for points in segment.points:
        currentElevation = points.elevation
        x.append(currentElevation)

        if currentElevation > prevElevation:
            elevgain+= (currentElevation - prevElevation)
            prevElevation = currentElevation
        else:
            prevElevation = currentElevation

    filt = np.convolve(x,kernel,mode='same')

    filt = np.delete(filt, 0)

    y = filt
    x = np.arange(len(filt))
    dydx = np.log10(abs(np.diff(y))+1)
    dydx = dydx - np.min(dydx)
    dydx = dydx / np.max(dydx)

    plt.plot(dydx)
    plt.show()

    fig, ax = plt.subplots()

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    norm = plt.Normalize(dydx.min(), dydx.max())
    lc = LineCollection(segments, cmap='plasma', norm=norm)

    lc.set_array(dydx)
    lc.set_linewidth(2)
    line = ax.add_collection(lc)
    fig.colorbar(line)

    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(y.min(), y.max())

    npts = len(y)
    colourmap = mpl.cm.get_cmap('plasma')
    normalize = mpl.colors.Normalize(vmin=dydx.min(), vmax=dydx.max())

    for i in range(npts - 1):
        plt.fill_between([x[i], x[i + 1]],
                         [y[i], y[i + 1]],
                         color=colourmap(normalize(dydx[i]))
                         , alpha=0.6)

    plt.show()

    for points in filt:
        currentFilteredElevation = points

        if currentFilteredElevation > prevFiltElevation:
            filtElevgain+= (currentFilteredElevation - prevFiltElevation)
            prevFiltElevation = currentFilteredElevation
        else:
            prevFiltElevation = currentFilteredElevation

    plt.plot(x)
    plt.plot(filt)
    plt.ylim(50,1100)
    plt.title(("Táv", round(n.length_3d(),1), "emelkedés", round(elevgain,0)))
    plt.savefig('kk_day'+str(fname)+'.png')

    plt.show()
    fname = fname + 1
    print("Táv", round(n.length_3d(),1), "emelkedés", round(elevgain,0)," filtered elevgain: ",filtElevgain)
    elevgain = 0

print("hello")