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
from collections import namedtuple, Counter
import cv2
import numpy as np
from scipy.spatial import ConvexHull, Delaunay
from threading import Thread, Condition, Lock
from Queue import Queue, Empty
import logging

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s (%(threadName)-2s) %(message)s',
                    )

# GLOBALS

DEBUG = True
counter = 0
colors = {}
tmpdir = "."
settings_stack = []
sessionfiles = []
images = []
image_mask = None
image_queue = Queue()
entropy_queue = Queue()
condition = Condition()
entropies = []
lock = Lock()

FALSE_COLOR_SETTINGS = {
#    'cartoon_discrete_colors': 1,
#    'cartoon_use_shader': 1,
    'orthoscopic': 1,
    'pick_shading': 1,
    'bg_rgb': 'black'
}


class ImageSampler:  

    def run(self, views, width, height, features, keepFile, callback = None):
        # init state
        self.views = views
        self.width = width
        self.height = height
        self.features = features
        self.keepFile = keepFile
        self.callback = callback
        self.pick_shading = 0
        self.task_id = 0
        self.results = []
       
        if len(self.views) == 0:
            print "need to set at least one view"
            return

        # setup callback
        cmd.raw_image_callback = self.capture_image

        # flat shading (patched PyMOL)
        self.pick_shading = cmd.get_setting_int('pick_shading')
        cmd.set('pick_shading')

        # start the recursion
        self.next()

    def cleanup(self):
        cmd.raw_image_callback = None
        cmd.set('pick_shading', False)
        
        if self.callback is not None:
            self.callback(self.features, self.results)

    def capture_image(self, img):
        '''
        callback, gets called after cmd.draw
        '''
        logging.debug('capturing view %i', self.task_id)
        e = self.image_entropy(img)
        self.results.append((self.views[self.task_id], e))

        self.task_id += 1
        self.next()
        
    def next(self):
        '''
        set up the next task and render one image
        '''

        if self.task_id >= len(self.views):
            self.cleanup()

            # report results
            logging.debug('done')
            return

        # set view
        logging.debug('setting view: %i', self.task_id)
        cmd.set_view(self.views[self.task_id])

        # render image
        cmd.draw(self.width, self.height, antialias=0)

        # self.task_id += 1


    def image_entropy(self, img):
        """compute the entropy of an image"""
        size = img.shape[0] * img.shape[1]

        img.resize((size, 4))
        img_flat = img[:,0] + img[:,1] * 256 + img[:,2] * 65536
        color_count = Counter(img_flat)
        

        if DEBUG:
            print "found %i different colors" % len(color_count)

        samples_probability = [float(count) / size for count in color_count.values()]

        e = self.entropy(samples_probability)

        if DEBUG:
            print "image entropy is: ", e

        return e
    
    def entropy(self, samples_probability):
        return -sum([p * log(p, 2) for p in samples_probability])



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
    feature_count = assign_colors(selection, by)

    apply_false_color_settings()
    compute_viewpoint_entropy(feature_count, view, width, height, keepFile)

    #pop()

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

    cmd.color('black', 'all')
    cmd.spectrum('b', 'rainbow', selection)

    if DEBUG:
        print "number of features: ", counter
    return counter

def print_results(features, results):

    for result in results:

        e = result[1]
        # normalize by number of features
        e /= log(features + 1, 2)
        print 'viewpoint entropy: ', e

# apply_false_color_settings() has to be called before running this function
def compute_viewpoint_entropy(features, view, width, height, keepFile=''):    
    logging.debug('compute_viewpoint_entropy')
    global ist

    if view:
        cmd.set_view(view)

    ist.run([view], width, height, features, False, print_results)

def capture_image(img):
    global entropy_queue
    global condition

    print 'image captured'
    return

    condition.acquire()
    e = image_entropy(img)
    print 'putting entropy ', e
    entropy_queue.put(e)
    condition.notify()
    condition.release()

def apply_false_color_settings():

    if DEBUG:
        print "setting configuration ", FALSE_COLOR_SETTINGS

    for key in FALSE_COLOR_SETTINGS:
        value = FALSE_COLOR_SETTINGS[key]
        cmd.set(key, value)

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

    feature_count = assign_colors(selection, by)
    apply_false_color_settings()

    keepFile = ''
    if DEBUG:
        keepFile = 'debug'

    #sampler = ViewpointEntropySampler(views, feature_count, width, height, keepFile)
    #sampler.sample()

    # return sampler.get_entropies()

    e = []
    for view in views:

        entropy = compute_viewpoint_entropy(feature_count, view, width, height, keepFile)
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


def best_view(selection='all', by='residues', n=10, width=100, height=100, ray=0, prefix='', add_PCA = False):
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
    features = assign_colors(selection, by)
    apply_false_color_settings()

    ist.run([view for (point, view) in views], width, height, features, False, set_best_view)


def set_best_view(features, results):
    '''
    callback to set the best view from image capture results
    '''
    ret = []

    maxi = -1.0;
    best_view = None
    for view, entropy in results:
       if entropy > maxi:
           maxi = entropy
           best_view = view

    cmd.set_view(best_view)
    pop()


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
        logging.debug("pushing %s", sessionfile)

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

        logging.debug("popping %s", sessionfile)
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
            cmd.set(key, settings[key])

        if DEBUG:
            print "restoring configuration ", settings


def cleanup():
    """remove all temporary data"""
    if not DEBUG:
        rmtree(tmpdir)

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

def PIL2array(img):
    return np.array(img.getdata(),
                    np.uint8).reshape(img.size[1], img.size[0], 3)

def array2PIL(arr, size):
    mode = 'RGBA'
    arr = arr.reshape(arr.shape[0]*arr.shape[1], arr.shape[2])
    if len(arr[0]) == 3:
        arr = np.c_[arr, 255*np.ones((len(arr),1), np.uint8)]
    return Image.frombuffer(mode, size, arr.tostring(), 'raw', mode, 0, 1)


ist = ImageSampler()

cmd.extend('best_view', best_view)

if DEBUG:
    cmd.extend('viewpoint_entropy', viewpoint_entropy)
    cmd.extend('assign_colors', assign_colors)
    cmd.extend('apply_false_color_settings', apply_false_color_settings)
    cmd.extend('push', push)
    cmd.extend('pop', pop)