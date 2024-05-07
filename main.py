import gpxpy
from gpxslicer import slicer
import matplotlib.pyplot as plt
import numpy as np
from gpxplotter import read_gpx_file, create_folium_map, add_segment_to_map
import folium
from math import pi, sqrt, exp


def gauss(n,sigma):
    r = range(-int(n/2),int(n/2)+1)
    return [1 / (sigma * sqrt(2*pi)) * exp(-float(x)**2/(2*sigma**2)) for x in r]

kernel = gauss(10,2)

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

days = 3

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

    webbrowser.open('KK_day'+str(fname)+'.html')

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

    filt = np.convolve(x,kernel,mode='valid')

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