from pymol import cmd, stored
from pymol.cgo import *
from PIL import Image
from tempfile import mkdtemp
from shutil import rmtree
from math import sin,cos,pi,sqrt,log,acos,degrees,acos,atan2
import os.path
import time
import string
from cgkit.cgtypes import *
import colorsys
import csv
from axes import *
import com #center of mass
import random
from collections import namedtuple
import cv2
import numpy as np
from scipy.spatial import ConvexHull, Delaunay

DEBUG = True

FALSE_COLOR_SETTINGS = {
    'specular': 0,
    'light_count': 1,
    'ambient': 1,
    'direct': 0,
    'shininess': 0,
    'reflect': 0,
    'depth_cue': 0,
    'fog': 0,
    'cartoon_discrete_colors': 1,
    'cartoon_use_shader': 1,
    'orthoscopic': 1
    # 'bg_rgb': 'black' # breaks with dangling pointer
}

def viewpoint_entropy(selection='all', by='residues', view=None, width=512, height=512, keepFile=''):
    '''
DESCRIPTION
 
    Computes the viewpoint entropy of 'selection' from the given 'view' by rendering the scene to a file of 
    the given dimensions. Keeps the file if 'keepFile' is set to true.
    If no 'view' is given, uses the current view.
    The viewpoint entropy can be computed either by secondary structure ('ss'), by 'residues' or by 'atoms'.


AUTHOR
 
    Julian Heinrich
    julian@joules.de
 
USAGE
 
    viewpoint_entropy selection=string, width=int, height=int
 
EXAMPLES
 
    viewpoint_entropy n.CA, 1024, 768
    '''
    
    push()

    # formalise parameters
    width, height = int(width), int(height)

    # get view
    if view is None:
        view = cmd.get_view(quiet=not DEBUG)

    # assign unique colors to each residue/atom for the selection
    assign_colors(selection, by)

    e = compute_viewpoint_entropy(view, width, height, keepFile)

    pop()

    if DEBUG:
        print "Viewpoint Entropy is: ", e
    return e


def assign_colors(selection, by='residues'):
    '''
    Assigns a unique color to every atom, residue, ss, or chain in the selection.
    '''
    stored.index = []
    stored.atoms = []
    myspace = {'atoms':[]}
    if by != 'atoms':
        selection += ' and n. CA'
    cmd.iterate(selection, 'stored.atoms.append((model,chain,ss,resi,name))')

    counter = 0
    increment = 1
    previous = 'hfjkdasnck12dd32' # pseudo-random

    if DEBUG:
        print "coloring by ", by

    for (model, chain, ss, resi, name) in stored.atoms:
        if by == 'atoms':
            counter += increment
        elif by == 'residues':
            if resi != previous:
                counter += increment
                previous = resi
        elif by == 'ss':
            if ss != previous:
                counter += increment
                previous = ss
        elif by == 'chain':
            if chain != previous:
                counter += increment
                previous = chain

        # alternating colors from both ends of the spectrum
        stored.index.append(counter if counter % 2 else -counter) 

    cmd.alter(selection, 'b = stored.index.pop(0)')
    cmd.spectrum('b', 'rainbow', selection)

    if DEBUG:
        print "number of features: ", counter
    return counter


def compute_viewpoint_entropy(features, view, width, height, keepFile=''):
    
    cmd.set_view(view)

#    apply_false_color_settings()
    tmpFile = render_false_color_image(width, height, keepFile)

    e = 0.0
    if os.path.isfile(tmpFile):
        img = Image.open(tmpFile)
        e = image_entropy(img)
        # normalize by number of features
        e /= log(features + 1, 2)

        if not len(keepFile):
            os.remove(tmpFile)

    return e



    

def apply_false_color_settings(useShader = 1):

    for key in FALSE_COLOR_SETTINGS:
        value = FALSE_COLOR_SETTINGS[key]
        cmd.set(key, value)

def render_false_color_image(width, height, filename=''):
    global tmpdir

    tmpFile = get_temporary_file(filename)

    # draw with antialiasing disabled
    cmd.draw(width, height, 0)
    cmd.png(tmpFile)
    while not os.path.exists(tmpFile):
        time.sleep(1)

    if DEBUG:
        print "created temporary file %s" % tmpFile

    return tmpFile

def get_temporary_file(filename=''):
    global tmpdir
    i = 0
    tmpFile = "%s/%s_%i.png" % (tmpdir, filename, i)
    while os.path.exists(tmpFile):
        i += 1
        tmpFile = "%s/%s_%i.png" % (tmpdir, filename, i)

    return tmpFile


def sample_viewpoint_entropy(views, selection='all', by='residues', width=512, height=512):

    # formalise parameters
    width, height = int(width), int(height)

    m = assign_colors(selection, by)

    keepFile = ''
    if DEBUG:
        keepFile = 'debug'

    e = []
    for view in views:

        entropy = compute_viewpoint_entropy(m, view, width, height, keepFile)
        e.append((view, entropy))
        
    return e

def get_views(points):
    ''' computes view matrices from points.
        assumes that the current view defines the up vector '''

    view = cmd.get_view(quiet=not DEBUG)
    rot = get_rotation(view)
    cam = rot.getRow(2).normalize()
    model_up = rot.getRow(1).normalize()
    model_right = rot.getRow(0).normalize()

    # temporarily decrease epsilon used for comparisons of vectors in cgit
    # to prevent math domain error when computing angles
    epsilon = getEpsilon()
    setEpsilon(9e-4)

    e = []
    for point in points:
        
        # find rotation matrix R to the new point
        spn = vec3(point).normalize()
        q = rotation_to(spn, cam)
        R = q.toMat3()

        # add rotation to initial rotation
        RR = rot * R

        # new camera coordinate system aligns
        # up vector with model_up vector
        # new_cam = vec3(point).normalize()
        new_cam = RR.getRow(2).normalize()
        lookAt = -new_cam
        up = model_right.cross(lookAt).normalize()

        #print "up: " + str(up) + "model_up: " + str(model_up)
        #print up == model_up
        #print up.cross(model_up).length() < 9e-4
        if up != model_up:
            if up.cross(model_up).length() < 9e-4 or up.angle(model_up) > pi/2:
                up = -(model_right).cross(lookAt).normalize()

        right = lookAt.cross(up).normalize()

        # compute the new view matrix
        RRR = mat3()
        RRR.setRow(0, right)
        RRR.setRow(1, up)
        RRR.setRow(2, new_cam)
        
        new_view = []
        new_view.extend(RRR.toList())
        new_view.extend(view[9:18])

        e.append((point, new_view))
        
    setEpsilon(epsilon)

    return e

def get_cam(view):
    ''' extracts the viewpoint (look-at) from the view matrix'''
    rot = get_rotation_from_view(view)
    cam = rot.getRow(2).normalize()

    return cam

def draw_up_vector(radius=1.0, rgb=[1.0, 1.0, 1.0], selection='all'):
    view = cmd.get_view(quiet=1)
    rot = get_rotation_from_view(view)
    cam = rot.getRow(2).normalize()
    model_up = rot.getRow(1).normalize()
    right = rot.getRow(0).normalize()

    bb = cmd.get_extent(quiet=not DEBUG, selection=selection)
    ll = vec3(bb[0])
    tr = vec3(bb[1])
    diag = tr - ll
    r = diag.length() / 4.0

    c = com.COM(selection)
    cvec = vec3(c)
    #cvec = vec3(view[12:15])

    drawVector(cvec, cvec + r * model_up, name='up')
    drawVector(cvec, cvec + r * right, name='right', rgb=[1.0, 0.0, 0.0])
    drawVector(cvec, cvec + r * cam, name='cam', rgb=[0.0, 1.0, 0.0])

def get_rotation(view):
    rot = mat3(view[0:9])
    # mat3 expects row-wise in constructor
    rot = rot.transpose()
    # now rot resembles view, i.e. column-major
    # matrix which rotates model axes to camera axes
    return rot


def best_view(selection='all', by='residues', n=10, width=512, height=512, ray=0, prefix='', add_PCA = False):
    # formalise parameters
    width, height, n = int(width), int(height), int(n)
    selection = str(selection)

    push()

    # sample points on sphere
    points = hammersley_points(n)
    
    if add_PCA:
        (pca1, pca2) = get_PCA_views(selection)
        points.append(pca1)
        points.append(pca2)

    # if DEBUG:
    #     # works only in this order for some reason...
    #     show_points(points, [1.0, 0.0, 0.0], selection)
    #     draw_up_vector(selection=selection)
    #     cmd.set("ray_trace_mode", '1')
    #     cmd.zoom('all', 2.0, 0, 1)

    views = get_views(points)

    # compute attributes
    entropies = sample_viewpoint_entropy([view for (point, view) in views], selection, by, width, height)

    # ret = []

    # basename = string.replace(selection," ","_")

    maxi = -1.0;
    best_view = None
    for view, entropy in entropies:
        if entropy > maxi:
            maxi = entropy
            best_view = view

    pop()    

    cmd.set_view(best_view)
    return best_view


def get_PCA_views(selection):
    old_view = cmd.get_view(quiet=1)

    cmd.orient(selection)
    view = cmd.get_view(quiet=1)
    rot = get_rotation_from_view(view)
    
    # now rot resembles view
    pc3 = rot.getRow(2)
    pc1 = rot.getRow(0)
    pc2 = rot.getRow(1)

    preferred1 = rot.getRow(2).normalize()

    rot2 = rot.rotate(pi, pc1)

    preferred2 = rot2.getRow(2).normalize()

    return (preferred1, preferred2)



def spherical_distance(a, b):
    avec = vec3(a).normalize()
    bvec = vec3(b).normalize()
    # return acos(avec.normalize() * bvec.normalize())
    return atan2(avec.cross(bvec).length(), avec * bvec)

def sample_and_show_points(selection='all', n=10, attribute = 'viewpoint_entropy', by='ss', width=512, height=512, keepFile = ''):
    #push()

    apply_false_color_settings()

    n = int(n)
    width, height = int(width), int(height)
    #points = []
    #(pca1, pca2) = getPCAViews(selection)
    #points.append(pca1)
    #points.append(pca2)
    points = hammersley_points(n) #poisson_disc(n, 50, points)

    views = get_views(points)
    attribs = []
    for (point, view) in views:
        a = viewpoint_attribute(attribute, selection, by, view, width, height, keepFile)
        attribs.append(a)

    show_points(points, [1.0, 0.0, 0.0], selection)
    show_attribute_sphere(points, attribs, selection, attribute)
    #pop()

def push_session():
    global tmpdir
    global sessionfiles

    if not os.path.isdir(tmpdir):
        if not DEBUG:
            tmpdir = mkdtemp()
        else:
            tmpdir = "."
 
    sessionfile = "%s/session_%i.pse" % (tmpdir, len(sessionfiles))

    if not os.path.isfile(sessionfile):
        cmd.save(sessionfile)
        sessionfiles.append(sessionfile)
        if DEBUG:
            print "pushing ", sessionfile

    # store current colors
    #stored.color_list = []
    # myspace = {'color_list':[]}
    #cmd.iterate(selection, 'stored.color_list.append((chain, resi, name, color))')

def pop_session():
    global sessionfiles

    if len(sessionfiles):
        sessionfile = sessionfiles.pop()
        if os.path.isfile(sessionfile):
            cmd.load(sessionfile)
            os.remove(sessionfile)
        else:
            print "SESSION STACK CORRUPTED!"

        if DEBUG:
            print "popping ", sessionfile
    else:
        cleanup()

    # also cleanup after final pop
    if not len(sessionfiles):
        cleanup()

def push():
    """save current state to disk"""
    global settings_stack

    push_session()

    settings = {}
    for key in FALSE_COLOR_SETTINGS:
        settings[key] = cmd.get(key)

    settings_stack.append(settings)

def pop():
    """restore last state that was saved via push"""
    global settings_stack

    pop_session()

    if len(settings_stack):
        settings = settings_stack.pop()
        for key in settings:
            print 'setting ', key, ' to ', settings[key]
            cmd.set(key, settings[key])

        if DEBUG:
            print "popping ", settings


def cleanup():
    """remove all temporary data"""
    if not DEBUG:
        rmtree(tmpdir)

def compute_image_entropy(image):
    global colors
    maxcolors = image.size[0]*image.size[1]
    print maxcolors
    imgColors = image.getcolors(maxcolors)
    print "found %i different colors for size %i" % (len(imgColors), maxcolors)

    used = 0
    e = 0
    N = log(maxcolors, 2)

    for count, color in imgColors:
        r = color[0]
        g = color[1]
        b = color[2]
        col_name = '0x%02x%02x%02x' % (r, g, b)
        if col_name in colors:
            if count > 0:
                # if DEBUG:
                #print "found " + str(count) + " pixels of colour " + col_name
                e = e + count * (log(count, 2) - N)
            used = used + 1
    
    e = -e / maxcolors
    #print 'e: ' , e, " from ", str(used), " features"
    return e

# assumes black background!!!
def image_area(image):
    size = image.size[0]*image.size[1]
    imgColors = image.getcolors(size)

    area = 0
    bg = (0, 0, 0, 255)
    for count, color in imgColors:
        if color != bg:
            area += count
    
    area = float(area)/float(size)
    ret = (area, len(imgColors))
    return ret

def image_entropy(image):
    """compute the entropy of an image"""

    size = image.size[0]*image.size[1]
    imgColors = image.getcolors(size)
    colorCount = len(imgColors)

    if DEBUG:
        print "found %i different colors" % colorCount

    samples_probability = [float(count) / size for (count, color) in imgColors]

    e = entropy(samples_probability)

    return e
    
def entropy(samples_probability):
    return -sum([p * log(p, 2) for p in samples_probability if p != 0])

def hammersley_points(n):
    '''computes hammersley points on the unit sphere'''
    points = []
    for k in range(n):
        t = 0
        p = 0.5
        kk = k
        while kk > 0:
            if (kk & 1):
                t += p
            p *= 0.5 
            kk >>= 1

        t = 2.0 * t - 1.0
        theta = (k + 0.5) / n
        thetarad = theta * 2.0 * pi     # theta in [0, 2pi]
        st = sqrt(1.0 - t*t)
        point = (st * cos(thetarad), st * sin(thetarad), t)
        #point = (thetarad, acos(t))
        points.append(point)

    return points;

# http://mathworld.wolfram.com/SphericalCoordinates.html
def to_spherical(x, y, z):
    r = sqrt(x*x + y*y + z*z)
    phirad = acos(z/r)
    thetarad = atan2(y,x)
    return (r, thetarad, phirad)

# http://mathworld.wolfram.com/SphericalCoordinates.html
def to_cartesian(r, thetarad, phirad):
    x = r * cos(thetarad) * sin(phirad)
    y = r * sin(thetarad) * sin(phirad)
    z = r * cos(phirad)
    return (x, y, z)

# values is assumed to be a list of point-value tuples
def show_points(points, color = [1.0, 0.0, 0.0], selection='all', name = 'samples', labels=[]):
    view = cmd.get_view(quiet=not DEBUG)

    # adjust radius to size of bounding box
    bb = cmd.get_extent(selection, quiet=not DEBUG)
    ll = vec3(bb[0])
    tr = vec3(bb[1])
    diag = tr - ll
    r = diag.length() / 2.0

    # origin of rotation in model space
    #o = vec3(view[12:15])
    c = com.COM(selection)
    o = vec3(c)
    #spheres = [BEGIN, TRIANGLE_STRIP]
    spheres = [COLOR]
    spheres.extend(color)
    i = 0.0
    j = 0
    for p in points:
        
        #spheres.extend([COLOR, 1.0, 1 - scaled_value, 1 - scaled_value])
        spheres.extend([SPHERE, o[0]+r*p[0], o[1]+r*p[1], o[2]+r*p[2], 1.25])
        #drawVector(o, o + r * vec3(p))
        #spheres.extend([VERTEX, o[0]+r*p[0], o[1]+r*p[1], o[2]+r*p[2]])
        #i += 1.0/len(values)
        l = 1.1
        if (len(labels) > j):
            cmd.pseudoatom(labels[j] + "_label", pos = [o[0]+l * r*p[0], o[1]+l*r*p[1], o[2]+l*r*p[2], 1.25], label = labels[j])    

        j += 1

    #spheres.extend([END])
    cmd.load_cgo(spheres, name, 1)

def get_rotation_from_view(view):
    rot = mat3(view[0:9])
    # mat3 expects row-wise in constructor
    rot = rot.transpose()
    # now rot resembles view, i.e. column-major
    # matrix which rotates model axes to camera axes
    return rot

def rotation_to(a, b):
    dot = a * b
    if dot < -0.99999999:
        tmpVec3 = vec3(1,0,0).cross(a)
        if tmpVec3.length() < 0.000001:
            tmpVec3 = vec3(0,1,0).cross(a)
        tmpVec3 = tmpVec3.normalize()
        out = quat(pi, tmpVec3)
        return out
    elif dot > 0.999999:
        return quat()
    else:
        tmpVec3 = a.cross(b)
        out = quat(1 + dot, tmpVec3[0], tmpVec3[1], tmpVec3[2])
        return out.normalize()

# FIXME: shouldn't use globals
counter = 0
colors = {}
tmpdir = "."
settings_stack = []
sessionfiles = []


cmd.extend('best_view', best_view)
if DEBUG:
    cmd.extend('assign_colors', assign_colors)
    cmd.extend('apply_false_color_settings', apply_false_color_settings)
    cmd.extend('push', push)
    cmd.extend('pop', pop)