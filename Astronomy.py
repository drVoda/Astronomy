import matplotlib.pyplot as plt
import glob
from math import *
from statistics import mean
import scipy.signal
from scipy import interpolate
from scipy.constants import G

def interpolation(x, x_points, y_points):
    tck = interpolate.splrep(x_points, y_points)
    return interpolate.splev(x, tck)

def smoothen(wave, n):
    return scipy.signal.savgol_filter(wave, n, 3)

def sort(s, n):
    for i in range(1, n):

        temp = s[i]
        j = i - 1

        while j >= 0 and len(temp) < len(s[j]):
            s[j + 1] = s[j]
            j -= 1

        s[j + 1] = temp

    return s

def bracket(start, end, size, x, y):
    x1 = []
    y1 = []
    n = int((end - start) / size)
    counter = start
    e = []
    for i in range(n):
        sx = 0
        sy = 0
        nx = 0
        ny = 0
        for j in range(len(x)):
            if x[j] > counter and x[j] < counter + size:
                sx += x[j]
                nx += 1
                sy += y[j]
                ny += 1
        x1.append(counter+size/2)
        y1.append(sy / ny)
        error = 0
        for j in range(len(x)):
            if x[j] > counter and x[j] < counter + size:
                if abs(y[j] - y1[i]) > error:
                    error = abs(y[j] - y1[i]) / y1[i]
        e.append(error)
        counter += size
    #print(e)
    return x1, y1

def data_analysis(path):
    files = glob.glob(path)
    files = sort(files, len(files))
    distance_from_galactic_center = []
    rotational_velocity = []
    vel_max = []
    longitude = []
    latitude = []

    for i in files:

        dat = open(i, 'r')
        vel = []
        temp = []
        c = False
        dataset_type = dat.readline()[0]
        dat.seek(0)

        for j in dat.readlines():

            if dataset_type == '#':

                if 'GLON=' in j:
                    longitude.append(radians(float(j[7:])))

                if 'GLAT=' in j:
                    latitude.append(radians(float(j[7:])))

                if j[0] != '#':
                    vel.append(float(j[:j.find(' ')]))
                    temp.append(float(j[j.find(' '):]))

            else:

                if '%%   ' in j:
                    longitude.append(radians(float(j[5:13])))
                    latitude.append(0)

                if '%%LAB' in j:
                    c = True

                if j[0] != '%' and c:
                    vel.append(float(j[:12]))
                    temp.append(float(j[12:]))

        temp_smoothened = smoothen(temp, 13)
        baseline = mean(temp_smoothened)
        v = -1000

        for j in range(1, len(temp_smoothened) - 2):

            if temp_smoothened[j] > baseline and temp_smoothened[j - 1] < temp_smoothened[j] and temp_smoothened[j + 1] < temp_smoothened[j] and vel[j] > v:
                v = vel[j]
                t = temp_smoothened[j]
                plt.plot(v, t, 'ro')

        vel_max.append(v)

        naslov='Spektar za l = 20° poslije primjenjivanja Savitzky - Golay algoritma'
        plt.plot(vel, temp_smoothened, 'b-', vel, [baseline]*len(vel), 'g--')
        plt.xlabel('Brzina relativna na LSR [Km / s]')
        plt.ylabel('Temperatura antene [K]')
        plt.title(naslov)
        plt.show()

    for i in range(len(longitude)):
        distance_from_galactic_center.append(8.5 * sin(longitude[i]))
        rotational_velocity.append(vel_max[i] / cos(latitude[i]) + 220 * sin(longitude[i]))

    d=distance_from_galactic_center[1:]
    r=rotational_velocity[1:]
    l=longitude[1:]
    #makeTable(['Longituda [rad]', 'Udaljenost [Kpc]', 'Brzina [Km / s]'], [l, d, r])

    return distance_from_galactic_center, rotational_velocity

def kepler(start):
    k = start
    i = 0.8
    x = [0]
    y = [0]
    while i < 8.5:
        x.append(i)
        y.append(k/sqrt(i))
        i += 0.01

    return x, y

def mass(radius, velocity):
    m = []
    for i in range (len(radius)):
        m.append(((radius[i]*3.08567758*(10**19))*((velocity[i]*1000)**2))/G)

    plt.plot(radius, m)
    plt.xlabel('Udaljenost od središta galaksije [kpc]')
    plt.ylabel('Masa galaksije [kg]')
    plt.title('Ovisnost mase o udaljenosti od središta')
    plt.show()
    return m

def makeTable(headerRow,columnizedData,columnSpacing=2):

    from numpy import array,max,vectorize

    cols = array(columnizedData,dtype=str)
    colSizes = [max(vectorize(len)(col)) for col in cols]

    header = ''
    rows = ['' for i in cols[0]]

    for i in range(0,len(headerRow)):
        if len(headerRow[i]) > colSizes[i]: colSizes[i]=len(headerRow[i])
        headerRow[i]+=' '*(colSizes[i]-len(headerRow[i]))
        header+=headerRow[i]
        if not i == len(headerRow)-1: header+=' '*columnSpacing

        for j in range(0,len(cols[i])):
            if len(cols[i][j]) < colSizes[i]:
                cols[i][j]+=' '*(colSizes[i]-len(cols[i][j])+columnSpacing)
            rows[j]+=cols[i][j]
            if not i == len(headerRow)-1: rows[j]+=' '*columnSpacing

    line = '-'*len(header)
    print(line)
    print(header)
    print(line)
    for row in rows: print(row)
    print(line)

def published():

    dat=open('C:\\Tekst3\objavljeno.txt')
    r=[]
    v=[]

    for i in dat.readlines():
        r.append(8.5*sin(radians(float(i.split(' ')[0]))))
        v.append(float(i.split(' ')[1])+220*sin(radians(float(i.split(' ')[0]))))

    #plt.plot(r, v, 'b.')
    #plt.show()
    return r, v

x, y = data_analysis('C:\\Data\*.txt')

x1 = []
y1 = []
i = 0
while i < 8.5:
    x1.append(i)
    y1.append(interpolation(i, x[1:], y[1:]))
    i += 0.01

x2, y2 = bracket(0, 1, 0.1, x1, y1)
x3, y3 = bracket(1, 8.5, 0.5, x1, y1)

x_RotationCurve = x2 + x3
y_RotationCurve = y2 + y3

x_kepler, y_kepler = kepler(interpolation(1, x[1:], y[1:]))

plt.plot(x, y, 'k.', x_RotationCurve, smoothen(y_RotationCurve, 5), 'r-') #, x_kepler, y_kepler, 'b-'
plt.xlabel('Udaljenost od središta galaksije [kpc]')
plt.ylabel('Obodna brzina [Km / s]')
plt.title('Rotacijski profil Mliječne staze')
plt.show()

mass(x_RotationCurve, y_RotationCurve)

#print(((8.5*3.08567758*(10**19))*((220*1000)**2))/G)

x_published, y_published = published()
x_publishedb, y_publishedb=bracket(2, 8.5, 0.5, x_published, y_published)
x_observedb, y_observedb = bracket(2, 8.5, 0.5, x1, y1)

for i in range(13):
    print(abs(y_publishedb[i]-y_observedb[i])/y_observedb[i], 2+i/2)
    print(y_observedb[i], y_publishedb[i], 2+i/2)

plt.plot(x_publishedb, y_publishedb, 'b.', x_observedb, y_observedb, 'k.')
plt.show()