from pymol import cmd, stored
from pymol.cgo import *
from tempfile import mkdtemp
from shutil import rmtree
from math import sin,cos,pi,sqrt,log,acos,degrees,acos,atan2,floor
import os.path
import time
import string
from cgtypes import *
import colorsys
from axes import *
import com #center of mass
import random
from collections import namedtuple, Counter
import numpy as np
import logging
from sklearn.cluster import MeanShift, estimate_bandwidth

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s (%(threadName)-2s) %(message)s',
                    )

# GLOBALS

DEBUG = True
tmpdir = "."
settings_stack = []
sessionfiles = []

# These settings are used for computing the viewpoint entropy
FALSE_COLOR_SETTINGS = {
    'surface_color_smoothing': 0,
    'orthoscopic': 1,
    'pick_shading': 1,
    'bg_rgb': 'black'
}

# Samples images from given views, width, and height.
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

        # flat shading (requires patched PyMOL)
        self.pick_shading = cmd.get_setting_int('pick_shading')
        cmd.set('pick_shading')

        # start recursion
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

        # if (len(self.keepFile)):
        #     im = Image.fromarray(img)
        #     im.save(self.keepFile + str(self.task_id) + ".png")

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
        """computes the entropy of an image"""
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

    Computes the viewpoint entropy of 'selection' from the given 'view' by rendering the scene to an image of
    the given dimensions.
    If no 'view' is given, uses the current view.
    The viewpoint entropy can be computed either by secondary structure ('ss'), by 'residues' or by 'atoms'.


AUTHOR

    Julian Heinrich
    julian@joules.de

USAGE

    viewpoint_entropy selection=string, by='residues|ss|atoms', view=[...], width=int, height=int

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
    visited = {}
    #if by != 'atoms':
    #    selection += ' and n. CA'
    cmd.iterate(selection, 'stored.atoms.append((model,chain,ss,resi,name))')

    counter = 0
    increment = 1
    previous = 'hfjkdasnck12dd32' # pseudo-random

    logging.debug("coloring by %s", by)

    for (model, chain, ss, resi, name) in stored.atoms:
        if by.startswith('atom'):
            counter += increment
        elif by.startswith('residue'):
            if resi != previous:
                counter += increment
                previous = resi
        elif by == 'ss':
            if ss != previous:
                counter += increment
                previous = ss
        elif by.startswith('chain'):
            if not chain in visited:
                counter += increment
                visited[chain] = True

        stored.index.append(counter)

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

    pop()

# apply_false_color_settings() has to be called before running this function
def compute_viewpoint_entropy(features, view, width, height, keepFile=''):
    global ist

    if view:
        cmd.set_view(view)

    ist.run([view], width, height, features, keepFile, print_results)


def apply_false_color_settings():
    '''
    applies FALSE_COLOR_SETTINGS
    '''
    logging.debug("setting configuration ", FALSE_COLOR_SETTINGS)

    for key in FALSE_COLOR_SETTINGS:
        value = FALSE_COLOR_SETTINGS[key]
        cmd.set(key, value)


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
        Assumes that the current view defines the up vector of the model '''

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

# Draws the up(white), right(red), and camera(green) vector of the current view.
def draw_up_vector(radius=1.0, rgb=[1.0, 1.0, 1.0], selection='all'):
    view = cmd.get_view()
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

def make_movie(views, n_frames = 100):
    cmd.mset("1 x100", 1)

    cmd.set_view(views[0])
    cmd.mview("store", 1)
    cmd.mview("store", 100)

    step = floor(n_frames/len(views))
    for i, view in enumerate(views[1:]):
        cmd.set_view(view)
        frame = (i + 1) * step
        logging.debug("setting view for frame " + str(frame))
        cmd.mview("store", frame)

    cmd.mplay()


def set_best_view_cb(features, results):
    '''
    callback to set the best view from image capture results
    also pops the settings stack
    '''
    ret = []

    (best_view, high_entropy) = max(results, key=lambda result: result[1])

    pop()
    cmd.set_view(best_view)

def show_scenes_cb(features, results):
    """
    Callback to create a set of scenes from the best views from image capture results.
    Also pops the settings stack.
    """

    # format entropies as array of arrays for estimate_bandwidth
    entropies = np.array([entropy for view, entropy in results]).reshape(-1, 1)
    bandwidth = estimate_bandwidth(entropies)
    ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
    ms.fit(entropies)
    labels = ms.labels_
    centers = ms.cluster_centers_.reshape(1, -1)[0]
    labels_unique = np.unique(labels)
    n_clusters = len(labels_unique)

    # find closest view to each center
    views = []
    for center in centers:
        d = 10 # maximum distance can be 1.0 by construction
        v = None
        for view, entropy in results:
            dd = abs(entropy - center)
            if dd < d:
                d = dd
                v = view
        views.append(v)

    # logging.debug(views)

    pop()

    cmd.set_view(views[0])
    cmd.scene('new', 'store')

    for i, view in enumerate(views[1:]):
        cmd.set_view(view)
        cmd.scene('new', 'store')


def loop_scenes():
    cmd.mset("1 x1", 1)
    cmd.set('scene_loop')
    cmd.set('movie_fps', 0.5)
    cmd.mdo(1, 'scene auto, next')
    cmd.mplay()


def show_orbit_cb(features, results):
    """Callback to orbit camera around object passing low- and high viewpoint entropy views.

    Keyword arguments:
    features -- not used
    results -- list of (view, entropy) pairs obtained from ImageSampler
    """

    # we want both the max and the min entropy to be on the camera orbit
    (view1, high_entropy) = max(results, key=lambda result: result[1])
    (view2, low_entropy) = min(results, key=lambda result: result[1])

    # find the vector orthogonal to the respective look at vectors
    look_at1 = get_cam(view1)
    look_at2 = get_cam(view2)

    # align z to the up-vector of the camera
    # compute the new view matrix
    up = look_at1.cross(look_at2).normalize()
    right = up.cross(look_at1).normalize()

    pop()

    R = mat3()
    R.setRow(0, right)
    R.setRow(1, up)
    R.setRow(2, look_at1)

    current_view = cmd.get_view()
    view = []
    view.extend(R.toList())
    view.extend(current_view[9:18])

    cmd.set_view(view)

    roll_camera(200, "y")


def roll_object (selection, frames, rotation_axis):
    angle = 360/frames
    frame = 1

    cmd.set("matrix_mode", 1)

    cmd.mset("1 x100", 1)

    cmd.frame(1)
    cmd.mview("store")
    cmd.mview("store", object=selection)

    while frame <= frames:
        cmd.frame(frame)
        cmd.rotate(rotation_axis, angle, object=selection)
        cmd.mview("store", frame)
        frame += 1

    cmd.mplay()

def roll_camera (frames, rotation_axis, filename = None):
    frames = int(frames)
    angle = 360.0/float(frames)
    frame = 1

    s = "1 x" + str(frames)
    cmd.mset(s, 1)

    cmd.mview("store", 1, power=1.0)

    while frame <= frames:
        cmd.turn(rotation_axis, angle)
        cmd.mview("store", frame, power = 1.0)
        if filename is not None:
            if ray:
                cmd.ray()
        frame += 1

    cmd.mplay()


def set_best_view(selection='all', by='residues', n=10, width=100, height=100, add_PCA = False, cb = set_best_view_cb):
    """
    Set the best view on a given selection.

    Keyword arguments:
    selection -- a PyMOL selection (default is 'all')
    by -- one of ('atoms', 'residues', 'ss', 'chain')
    n -- the number of samples to evaluate
    width -- the width (in pixels) for each sample
    height -- the height (in pixels) for each sample
    add_PCA -- include (TRUE) PCA viewpoints in addition to the n samples? (default is False)
    cb -- a callback reserved for internal use
    """

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

    views = get_views(points)
    features = assign_colors(selection, by)
    apply_false_color_settings()

    # TODO: call 'deselect' to remove all selections prior to rendering
    # also call 'zoom'

    # run image sampler for all viewpoints with cb as callback
    #cb = transition_to_best_view_cb if animate else set_best_view_cb
    ist.run([view for (point, view) in views], width, height, features, False, cb)

def show_scenes(selection='all', by='residues', n=10, width=100, height=100):
    """
    Computes a set of scenes from clustering viewpoint entropy.

    Keyword arguments:
    selection -- a PyMOL selection (default is 'all')
    by -- one of ('atoms', 'residues', 'ss', 'chain')
    n -- the number of samples to evaluate
    width -- the width (in pixels) for each sample
    height -- the height (in pixels) for each sample
    """
    set_best_view(selection, by, n, width, height, False, cb=show_scenes_cb)

def play_tour(features, results):
    """
    Callback to loop through best views.
    """
    show_scenes_cb(features, results)
    loop_scenes()

def tour(selection='all', by='residues', n=10, width=100, height=100):
    """
    Plays a camera tour passing scenes from clustering viewpoint entropy.

    Keyword arguments:
    selection -- a PyMOL selection (default is 'all')
    by -- one of ('atoms', 'residues', 'ss', 'chain')
    n -- the number of samples to evaluate
    width -- the width (in pixels) for each sample
    height -- the height (in pixels) for each sample
    """
    set_best_view(selection, by, n, width, height, False, cb=play_tour)

def orbit(selection='all', by='residues', n=10, width=100, height=100):
    """
    Shows a 360-degree rotation around the selection passing the highest and lowest viewpoint-entropy views.

    Keyword arguments:
    selection -- a PyMOL selection (default is 'all')
    by -- one of ('atoms', 'residues', 'ss', 'chain')
    n -- the number of samples to evaluate
    width -- the width (in pixels) for each sample
    height -- the height (in pixels) for each sample
    """
    set_best_view(selection, by, n, width, height, False, cb=show_orbit_cb)

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

    return atan2(avec.cross(bvec).length(), avec * bvec)

def push_session():
    '''
    pushes the current session on a stack and saves it to a file
    '''
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

def pop_session():
    '''
    pops the session from the stack and removes the corresponding file
    '''
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

        logging.debug("restoring configuration")


def cleanup():
    """removes all temporary data"""
    if not DEBUG:
        rmtree(tmpdir)

def hammersley_points(n):
    '''computes n hammersley points on the unit sphere'''
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
    '''converts Cartesian to spherical coordinates'''
    r = sqrt(x*x + y*y + z*z)
    phirad = acos(z/r)
    thetarad = atan2(y,x)
    return (r, thetarad, phirad)

# http://mathworld.wolfram.com/SphericalCoordinates.html
def to_cartesian(r, thetarad, phirad):
    '''converts spherical to Cartesian coordinates'''
    x = r * cos(thetarad) * sin(phirad)
    y = r * sin(thetarad) * sin(phirad)
    z = r * cos(phirad)
    return (x, y, z)

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

ist = ImageSampler()

cmd.extend('set_best_view', set_best_view)
cmd.extend('show_scenes', show_scenes)
cmd.extend('tour', tour)
cmd.extend('orbit', orbit)
cmd.extend('loop_scenes', loop_scenes)

if DEBUG:
    cmd.extend('viewpoint_entropy', viewpoint_entropy)
    cmd.extend('assign_colors', assign_colors)
    cmd.extend('apply_false_color_settings', apply_false_color_settings)
    cmd.extend('push', push)
    cmd.extend('pop', pop)
    cmd.extend('draw_up_vector', draw_up_vector)
    cmd.extend('roll_camera', roll_camera)
