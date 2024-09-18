import subprocess
from pathlib import Path

import mrcfile
import numpy as np
import scipy
from skimage.measure import block_reduce
from skimage.morphology import ball, remove_small_objects
from skimage.transform import resize


def extract(image, boxes, boxsize):
    """
    Crop particles from provided image.

    Takes ZYX-image and ZYX coordinates (boxes), and 1D boxsize in px as input.

    Returns 4D stack (n, z, y, x)

    """
    particles = np.zeros((len(boxes), boxsize, boxsize, boxsize))

    i = 0

    for box in boxes:

        z, y, x = box[0], box[1], box[2]

        # extract from ZYX order
        particles[i, :, :, :] = image[int(z-boxsize/2):int(z+boxsize/2),
                                      int(y-boxsize/2):int(y+boxsize/2),
                                      int(x-boxsize/2):int(x+boxsize/2)]

        i = i+1

    return particles


class Picker:
    """Picks blobs.

    Inputs:
    tomogram: Path of reconstruction
    output: Path of output folder
    radius: Particle radius (A)
    override_angpix: specify A/pix instead of using information from header
    contamination_binning: bin to this level for contamination mask

    Functions:
    - getcont() returns contamination mask
    - minima_extract() returns coordinates of local minima
    - detect() returns cleaned coordinates and cropped particles from filtered tomogram
    - prefilt() returns statistics about raw pick signal
    - filt() filters picks based on signal statistics at these coordinates

    """

    def __init__(self,
                 tomogram: Path,
                 output: Path,
                 radius: int = 100,
                 contamination_binning: int = 1,
                 override_angpix: float = None): #noqa: RUF013

        self.name = tomogram.stem
        self.output = output
        self.radius = radius
        self.contamination_binning = contamination_binning

        self.rec = mrcfile.read(tomogram)

        if override_angpix is None:
            with mrcfile.mmap(tomogram, mode='r') as f:
                self.pixelsize = float(f.voxel_size.x)
        else:
            self.pixelsize = override_angpix

        # box size for checking whether pick contains density
        self.tilesize = int(3 * radius / self.pixelsize)
        if self.tilesize % 2 > 0:
            self.tilesize += 1

        # HP filter with strict cutoff
        radius1 = 0.001
        sigma1 = 0

        # LP filter to 0.5x particle radius
        radius2 = 0.5 * self.pixelsize / self.radius
        sigma2 = 0.001

        # Create lowpass-filtered version of the tomogram
        # TODO: do this also in python to remove dependencies
        # https://scikit-image.org/docs/stable/auto_examples/filters/plot_dog.html

        subprocess.run(['mtffilter', '-3dfilter',
                        '-radius1', str(radius1),
                        '-hi', str(sigma1),
                        '-l', f'{radius2},{sigma2}',
                        tomogram,
                        output / f'{self.name}_bp.mrc'],
                       stdout=subprocess.DEVNULL)

        self.lowres = mrcfile.read(output / f'{self.name}_bp.mrc')

    def getcont(self,
                gaussian=True,
                sigma=200,
                stdtimes=3,
                min_size=125,
                dilation=100):
        """Get contamination mask.

        Inputs:
        gaussian: remove low-frequency variations before thresholding?
        sigma: sigma for gaussian background filter (A)
        stdtimes: how many SD below mean are considered contamination?
        min_size: only keep contaminations of this size or bigger (px)
        dilation: dilate contaminatins by this much (A)

        Returns
        -------
        area: mask for contaminations

        """
        # First, do binning for contamination mask if requested

        if self.contamination_binning != 1:
            G = block_reduce(self.rec,
                             block_size = self.contamination_binning,
                             func=np.mean,
                             cval=self.rec.mean())

        else:
            G = self.rec

        # Gaussian filter to remove background from image
        if gaussian:

            # If tomogram was converted to byte, values might be saturated
            # Thus, convert to float32
            G = G.astype(np.float32)

            sigma_px = sigma/(self.pixelsize*self.contamination_binning)

            # Create volume capturing low-frequency background changes
            Gf = scipy.ndimage.gaussian_filter(G, sigma=sigma_px)

            # Subtract from original image

            G = G-Gf

        maskthres=G.mean()-stdtimes*G.std()

        mask = G < maskthres

        cmask = scipy.ndimage.binary_opening(mask, ball(1))

        clean = remove_small_objects(cmask, min_size=min_size)

        segmentation = scipy.ndimage.binary_closing(clean, ball(2))

        contamination_dilation = int(dilation/
                                     self.contamination_binning/
                                     self.pixelsize)

        area = scipy.ndimage.binary_dilation(segmentation, ball(contamination_dilation))

        area_rescaled = resize(area,
                               self.rec.shape,
                               cval=False)

        mrcfile.write(self.output / f'{self.name}_contamination.mrc',
                      area_rescaled.astype(np.int16),
                      overwrite=True)

        return area_rescaled

    def minima_extract(self,radius_times=4,inhibit=False):
        """Extract minima.

        Inputs:
        radius_times: how many radii should there be between picks?
        inhibit: pass inhibit=True to get more/closer picks

        Returns
        -------
        points: coordinates of minima, in order of array (ZYX)

        """
        locality = int(radius_times*self.radius / self.pixelsize)

        # If inhibit is False, just return all minima
        if not inhibit:
            minima = (
                self.lowres == scipy.ndimage.minimum_filter(self.lowres, locality)
            ).nonzero()
            points=[]
            for i in range(minima[0].shape[0]):
                points.append([minima[0][i],minima[1][i],minima[2][i]])
            return points

        # Otherwise, make a bit fancier minima extraction
        # Algorithm explained in S1 of publication

        locality2 = int(radius_times*self.radius/2 / self.pixelsize)

        minima = (
            self.lowres == scipy.ndimage.minimum_filter(self.lowres, locality2)
        ).nonzero()

        rawpoints={}

        for i in range(minima[0].shape[0]):
            z,y,x=minima[0][i],minima[1][i],minima[2][i]
            v=self.lowres[z,y,x]
            if v not in rawpoints:
                rawpoints[v]=[]
            rawpoints[v].append([[z,y,x],False])

        editmax=self.lowres.max()

        edited=self.lowres+(1-
                            (self.lowres == scipy.ndimage.minimum_filter(
                                self.lowres,locality2))
                            )*(editmax-self.lowres)

        inhibited=1

        while inhibited>0:

            inhibited=0
            minf=scipy.ndimage.minimum_filter(edited, locality)

            editmax=edited.max()

            for k in rawpoints.keys():
                for i in range(len(rawpoints[k])):
                    p=rawpoints[k][i]
                    minimum=minf[p[0][0],p[0][1],p[0][2]]
                    if minimum==k:
                        p[1]=False
                    elif minimum in rawpoints:
                        p[1]=minimum
                        inhibited+=1
                    else:
                        p[1]=True
                        inhibited+=1
                        z,y,x=p[0][0]-int(locality/2),p[0][1]-int(locality/2),p[0][2]-int(locality/2)
                        a,b,c=z+locality,y+locality,x+locality
                        x=0 if x<0 else x
                        y=0 if y<0 else y
                        z=0 if z<0 else z
                        a=edited.shape[0] if a>edited.shape[0] else a
                        b=edited.shape[1] if b>edited.shape[1] else b
                        c=edited.shape[2] if c>edited.shape[2] else c
                        edited[z:a,y:b,x:c]+=(edited[z:a,y:b,x:c]==minimum)*(editmax-edited[z:a,y:b,x:c])

            deleting={}
            for k in rawpoints.keys():
                for i in range(len(rawpoints[k])):
                    p=rawpoints[k][i]
                    if not isinstance(p[1], bool):
                        d=False
                        for j in rawpoints[p[1]]:
                            if isinstance(j[1], bool) and not j[1]:
                                d=True
                                break
                        if d:
                            if k not in deleting:
                                deleting[k]=[]
                            deleting[k].append(i)


            for k in deleting.keys():
                for i in range(len(deleting[k])-1,-1,-1):
                    point=rawpoints[k][deleting[k][i]][0]
                    edited[point[0],point[1],point[2]]=editmax
                    rawpoints[k].pop(deleting[k][i])
                if len(rawpoints[k])==0:
                    del rawpoints[k]
        points=[]
        for k in rawpoints.keys():
            for p in rawpoints[k]:
                points.append(p[0])
        return points


    def detect(self,area,radius_times=4,inhibit=False,detection_width=128):
        """Detect particles.

        Inputs:
        area: contamination-mask
        radius_times: how close can picks be (multiples of radius)
        inhibit: pass True to get more picks by using iterative strategy
        detection_width: only consider picks this many pixels from center in Z (px)

        Returns
        -------
        boxes: list of particle coordinates in ZYX order
        raw_particles: images of particles with box tile_sizes in N*ZYX order

        """
        # Get minima
        points=self.minima_extract(radius_times=radius_times,inhibit=inhibit)

        # Turn to coordiantes
        boxes = []
        for i in range(len(points)):

            # Current point coordinate
            z = points[i][0]
            y = points[i][1]
            x = points[i][2]

            # Check overlap with contamination map
            clean = not area[z, y, x]

            # Check whether at edge
            inside = (
                x - self.tilesize / 2 >= 0
                and x < self.rec.shape[2] - self.tilesize/2 + 1
                and y- self.tilesize / 2 >= 0
                and y < self.rec.shape[1] - self.tilesize/2 + 1
                and z- self.tilesize / 2 >= 0
                and z < self.rec.shape[0] - self.tilesize/2 + 1

                # in the original implementation, this is y / shape[1]
                and z <= self.rec.shape[0]*0.5+detection_width
                and z >= self.rec.shape[0]*0.5-detection_width
            )

            if clean and inside:
                boxes.append([z, y, x])

        raw_particles = extract(self.lowres, boxes, self.tilesize)

        return boxes, raw_particles


    def prefilt(self, raw_particles, stdtimes=1):
        """Calculate parameters for particle filtering based on signal.

        Inputs:
        raw_particles: particle images in ZYX order
        stdtimes: how many SD

        Returns
        -------
        particles_metrics: SD of foreground, SD of background for each particle
        stdthreshold: absolute threshold for foreground SD based on stdtimes

        """
        # Define Foreground as central circle in box
        x, y= np.mgrid[0:self.tilesize, 0:self.tilesize] - self.tilesize / 2 + 0.5
        condition2d = np.sqrt(x*x+y*y) > self.radius / self.pixelsize

        # z-project particles
        raw_particles2d=raw_particles.sum(axis=1)

        # Calculate foreground vs. background SD
        particles_metrics = np.zeros([raw_particles.shape[0],2])

        for p in range(raw_particles.shape[0]):
            raw = raw_particles2d[p, :, :]
            background = np.extract(condition2d, raw)
            foreground = np.extract(np.logical_not(condition2d), raw)
            particles_metrics[p,0] =foreground.std()
            particles_metrics[p,1] =background.std()

        # Foreground SD mean and SD
        stdmean,stdstd=particles_metrics[:,0].mean(),particles_metrics[:,0].std()

        # Define threshold for foreground SD
        stdthreshold=stdmean+stdstd*stdtimes

        return particles_metrics,stdthreshold

    def filt(self,boxes,particles_metrics,stdthreshold,remove_edge):
        """
        Filter particle picks based on image characteristics.

        Inputs:
        boxes: list of pick coordinates, NZYX-ordered.
        particles_metrics: Foreground/background SD from Picker.prefilt()
        stdthreshold: actual SD threshold, from Picker.prefilt()
        remove_edge: remove particles with higher background than foreground SD?

        Returns
        -------
        boxs_XYZ: XYZ coordiantes of retained particles.

        Also writes out coordinate file in XYZ order to output folder.

        """
        boxs_ZYX = []

        rawboxs={}

        # Keep only particles which have a foreground SD higher than stdthreshold
        # and higher foreground SD than background SD

        for p in range(particles_metrics.shape[0]):
            std=particles_metrics[p,0]
            if std >= stdthreshold and (std>particles_metrics[p,1] or
                                        (not remove_edge)):
                if std not in rawboxs:
                    rawboxs[std]=[]
                rawboxs[std].append([boxes[p][0],boxes[p][1],boxes[p][2]])

        rawboxs=[rawboxs[k] for k in sorted(rawboxs.keys())]
        for b in rawboxs:
            boxs_ZYX+=b

        # Transform from ZYX to XYZ
        boxs_ZYX =  np.stack(boxs_ZYX)
        boxs_XYZ = boxs_ZYX[:,[2,1,0]]

        # Write out coordinates as txt file
        with open(self.output/f"{self.name}.coords", "w") as f:
            np.savetxt(f, boxs_XYZ, delimiter=' ', fmt='%s %s %s')

        return boxs_XYZ
