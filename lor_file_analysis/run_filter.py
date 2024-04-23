import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib_scalebar.scalebar import ScaleBar
from scipy import optimize


def fix_npy_ending(file_loc):
    # file_loc: str - path of the file in question

    # adds '.npy' if the string does not end in '.npy'
    if not file_loc[-4:] == ".npy":
        file_loc += ".npy"
    return file_loc


def keplerslicer(volume, sliceplane, slicevalue):
    # first step is to define a new 2d plane that can be graphed
    # if the slice plane is "x" then the x values are held
    # constant and a new array of y ans z is created
    # Symmetric for all axis. x->yz, y->xz, z->xy
    if sliceplane == "x" or sliceplane == "X":
        useslice = volume[slicevalue]
    elif sliceplane == "y" or sliceplane == "Y":
        useslice = volume[:, slicevalue]
    elif sliceplane == "z" or sliceplane == "Z":
        useslice = volume[:, :, slicevalue]
    else:
        print("only 'x', 'y', 'z' allowed for `sliceplace`. Try again.")
        exit()
    return useslice


# shows what the zernike error from the image is
def zernikeerror(p, image, rows, cols, order):
    # cart = RZern(order); #something from the Zernike polynoial library
    cart = None

    # defines an object? Polynomial up to order?
    x = p[0]
    y = p[1]
    r = p[2]

    # makes a set of positions from left to right and top to bottom to check
    # same thing done for y
    ddx = np.linspace(-x / r, (rows - x) / r, rows)
    # makes a range of numbers of size rows from left to right
    ddy = np.linspace(-y / r, (cols - y) / r, cols)
    xv, yv = np.meshgrid(ddx, ddy)  # makes an x and y matrix of values
    cart.make_cart_grid(xv, yv, unit_circle=False)
    eval = cart.eval_grid(cart.fit_cart_grid(image)[0], matrix=True)
    # evaluation of Zernike Polynomials somehow
    # print(np.size(eval))
    return np.reshape(image - eval, rows * cols)


# takes an np array or location of np array and fits a zernike up to `order` then subtracts to create a filtered image
def plot_zernike_fitted(
    file, order, sliceplane, slicevalue, display=True, padding=0, padding_value=0
):
    # file: str - location of .npy file, np.array - np array object
    # order: int - highest order of Zernike function to fit
    # sliceplane:
    # slicevalue:
    # display: bool - displays all images as default
    # padding: int - number of padding layers to add to image file. helps with corner behavior of the zernike fit

    # uses `file` data type to either load as np array or to treat as np array
    if isinstance(file, (str,)):
        currentimage = np.load(fix_npy_ending(file))
    elif isinstance(file, (np.ndarray,)):
        currentimage = file
    else:
        print('File type " %s " not supported' % type(file))
        exit()

    # checks to see if needs to call keplerslicer to render 3d image into 2d array
    if len(np.shape(file)) == 3:
        img = keplerslicer(currentimage, sliceplane, slicevalue)
    elif len(np.shape(file)) == 2:
        img = currentimage
    else:
        print(
            "Dimensions of %s not allowed. Try again with 2D or 3D numpy array." % file
        )
        exit()

    img = np.pad(img, pad_width=padding, mode="constant", constant_values=padding_value)

    # assigns rows, cols the lengths in x, y of img
    rows, cols = img.shape

    # defines the radius to create unitcircle
    r0 = cols / 2.0 if rows > cols else rows / 2.0

    x0 = rows / 2.0  # x0 = half width
    y0 = cols / 2.0  # y0 = half height

    # initial vector
    p0 = [x0, y0, r0]

    # finds zernikes up to `order` that best fit img
    p1, success = optimize.leastsq(zernikeerror, p0, args=(img, rows, cols, order))
    # cart = RZern(order);
    cart = None

    # assigns final vector components to x, y, r
    x, y, r = p1[0], p1[1], p1[2]
    print("Best estimate coordinate origin: (x,y,r)=(%s, %s, %s)" % (x, y, r))

    ddx = np.linspace(-x / r, (rows - x) / r, rows)
    ddy = np.linspace(-y / r, (cols - y) / r, cols)

    xv, yv = np.meshgrid(ddx, ddy)
    cart.make_cart_grid(xv, yv, unit_circle=False)

    # creates numpy arrays of all images
    eval = cart.eval_grid(cart.fit_cart_grid(img)[0], matrix=True)
    fitted_image = (img - eval)[padding : rows - padding, padding : cols - padding]
    original_image = img

    # [padding:rows-padding, padding:cols-padding]

    if display:
        # plot zernike
        plt.figure(1)
        aplot = plt.imshow(eval, cmap="gray")
        plt.colorbar(aplot, orientation="vertical")

        # plot fitted graph
        plt.figure(2)
        bplot = plt.imshow(fitted_image, cmap="gray")
        plt.colorbar(bplot, orientation="vertical")

        # plot original graph
        # print(original_image[81][125])
        # print(original_image[109][109])
        plt.figure(3)
        cplot = plt.imshow(original_image, cmap="gray")
        plt.colorbar(cplot, orientation="vertical")
        # plt.grid(color='w')
        plt.show()

    return fitted_image


# simply displays a 3d or 2d .npy file as a 2d image
def show_image(source_file, pixels, slice=10, super_slice=0, display=True, vmax=None):
    # source_file: str - location of .npy file
    # pixels: int - number of pixels of one side of square np array
    # slice: int -
    # super_slice: int
    # display: bool - True runs plt.show() to display image
    # vmax: float - sets the maximum scale of the values in each pixel, effectively changing the contrast

    source = np.load(fix_npy_ending(source_file))

    print("Maximum pixel value: " + str(np.amax(source)))

    print("Dimensions of source: %s" % str(np.shape(source)))

    # tests if source is a 2d or 3d image
    if len(np.shape(source)) == 3:
        image = source[:, :, slice]
        if super_slice > 0:
            for i in range(super_slice):
                image += source[:, :, slice + i - int(super_slice / 2)]
    elif len(np.shape(source)) == 2:
        image = source.reshape((pixels, pixels))
    else:
        print(
            "Dimensions of %s not allowed. Try again with 2D or 3D numpy array."
            % source_file
        )
        exit()

    # displays source
    if display:
        fig, ax = plt.subplots()
        ax_img = ax.imshow(image, cmap="gray", vmax=vmax, extent=[-10, 10, -10, 10])
        ax.set_xlabel("Position (cm)", fontdict={"size": 12})
        ax.set_ylabel("Position (cm)", fontdict={"size": 12})
        color_bar = plt.colorbar(ax_img)
        scalebar = ScaleBar(
            1,
            units="cm",
            dimension="si-length",
            # label=None,
            # length_fraction=.05,
            # height_fraction=0.004,
            # width_fraction=None,
            # location=None,
            pad=0.5,
            border_pad=1.25,
            # sep=None,
            # frameon=None,
            # color=None,
            # box_color=None,
            # box_alpha=None,
            # scale_loc=None,
            # label_loc=None,
            # font_properties=None,
            # label_formatter=None,
            # scale_formatter=None,
            # fixed_value=None,
            # fixed_units=None,
            # animated=False,
            # rotation=None,
        )
        ax.add_artist(scalebar)
        plt.show()

    return image


def fft(image, filterloc, display=True):
    # image: np.array - np array of the image to be filtered
    # filterloc: str - path to the np array of the filter
    # display: bool - displays image as default

    image_F = np.fft.fftn(image)

    filter = np.load(fix_npy_ending(filterloc))

    symmetry = 2

    filter = 1 * filter

    cutoff_low = symmetry
    cutoff_high = pixels - symmetry
    # filter[cutoff_low:cutoff_high,cutoff_low:cutoff_high] = 1
    # filter[cutoff_high:,:] = 1
    # filter[:,cutoff_high:] = 1
    # filter[:cutoff_low, :] = 1
    # filter[:, :cutoff_low] = 1

    filter[cutoff_low:cutoff_high, :] = 1
    filter[:, cutoff_low:cutoff_high] = 1

    # now filter out the edge behavior

    fig, ax = plt.subplots()
    # fig, ay = plt.subplots()
    ax.imshow(np.real(filter))
    # ay.imshow(np.real(image_F))
    # plt.show()
    print("Dimensions of filter: %s" % str(filter.shape))
    cleaned_F = np.multiply(image_F, filter)

    symmetry = 30

    cutoff_low = symmetry
    cutoff_high = pixels - symmetry

    cleaned_F[cutoff_low:cutoff_high, :] = 1
    cleaned_F[:, cutoff_low:cutoff_high] = 1

    cleaned = np.fft.ifftn(cleaned_F)

    fig, az = plt.subplots()

    cleaned = np.real(cleaned)
    # cleaned = np.where(cleaned > 0,cleaned, 0.0)

    # conditionally displays image
    if display:
        az_ = az.imshow(cleaned, cmap="gray")
        az_bar = plt.colorbar(az_)
        plt.show()

    # np.save(source_file + "_clean", cleaned)
    return cleaned


def signal_to_noise(image, patch_tuple, display=True):
    patch_lowerleft, patch_upperright = patch_tuple

    patch_width, patch_height = np.subtract(
        np.asarray(patch_upperright), np.asarray(patch_lowerleft)
    )

    patch_average = np.average(image[patch_lowerleft, patch_upperright])
    patch_sigma = np.std(image[patch_lowerleft, patch_upperright])
    snr = patch_average / patch_sigma

    if display:
        fig, ax = plt.subplots()
        aplot = ax.imshow(image, cmap="gray")
        fig.colorbar(aplot, orientation="vertical")
        patch_rect = patches.Rectangle(
            patch_lowerleft,
            patch_width,
            patch_height,
            linewidth=1,
            edgecolor="r",
            facecolor="none",
        )
        ax.add_patch(patch_rect)
        plt.show()

        print(
            f"""--------------------------------------
signal-to-noise: {snr}
mean: {patch_average}
standard deviation: {patch_sigma}
--------------------------------------"""
        )

    return snr, patch_average, patch_sigma


def contrast_to_noise(image, background_patch, voi_patch, display=True, voi_avg=None):
    # get info for the background
    background_snr, background_avg, background_sigma = signal_to_noise(
        image, background_patch, display=False
    )

    if voi_avg is None:
        # get info for the voi
        voi_snr, voi_avg, voi_sigma = signal_to_noise(image, voi_patch, display=False)

    if display:
        background_lowerleft, background_upperright = background_patch
        voi_lowerleft, voi_upperright = voi_patch

        background_width, background_height = np.subtract(
            np.asarray(background_upperright), np.asarray(background_lowerleft)
        )
        voi_width, voi_height = np.subtract(
            np.asarray(voi_upperright), np.asarray(voi_lowerleft)
        )

        fig, ax = plt.subplots()
        aplot = ax.imshow(image, cmap="gray")
        fig.colorbar(aplot, orientation="vertical")
        background_rect = patches.Rectangle(
            background_lowerleft,
            background_width,
            background_height,
            linewidth=1,
            edgecolor="r",
            facecolor="none",
        )
        voi_rect = patches.Rectangle(
            voi_lowerleft,
            voi_width,
            voi_height,
            linewidth=1,
            edgecolor="g",
            facecolor="none",
        )
        ax.add_patch(background_rect)
        ax.add_patch(voi_rect)
        plt.show()

    cnr = (voi_avg - background_avg) / background_sigma
    print(f"contrast-to-noise ratio: {cnr}")
    return cnr


def get_tumor_region_non_functional(img, center, diameter, z_plane, total_dimensions):
    # center is defined in xcat software as the slices in x, y, and z
    # z_plane is defined in amide as the mm distance from the center of the phantom to the z-plane

    rows, cols = img.shape

    radius = diameter / 2
    mm_to_xcat_ratio = 1.5
    total_dimensions_mm = np.array(total_dimensions) * mm_to_xcat_ratio
    center = mm_to_xcat_ratio * np.array(center)
    center_corrected = center - 0.5 * total_dimensions_mm
    print(center_corrected)
    height = np.abs(center_corrected[2] - z_plane)
    new_radius = np.sqrt(radius * radius - height * height)

    image_to_mm_ratio = 2.0
    new_radius = image_to_mm_ratio * new_radius

    center_2d = image_to_mm_ratio * center_corrected[0:2] + np.array(
        [rows / 2, cols / 2]
    )
    print(tuple(center_2d))
    print(new_radius)

    def circle(x, y):
        f = (x - center_2d[0]) ** 2 + (y - center_2d[1]) ** 2 - new_radius**2
        return f

    total_sum, total_count, i, j = 0, 0, 0, 0
    rows, cols = img.shape

    while i < rows:
        while j < cols:
            if circle(i, j) <= 0:
                total_sum += img[i, j]
                total_count += 1
            j += 1
        i += 1
    avg = total_sum / total_count

    fig, ax = plt.subplots()
    aplot = ax.imshow(img, cmap="gray")
    fig.colorbar(aplot, orientation="vertical")
    voi_circle = patches.Circle(
        xy=center_2d, radius=new_radius, linewidth=1, edgecolor="g", facecolor="none"
    )
    ax.add_patch(voi_circle)
    plt.show()

    return avg


def circle(x, y, center, radius):
    f = (x - center[0]) ** 2 + (y - center[1]) ** 2 - radius**2
    # print(f)
    return f


def get_tumor_region(img, center, edge, display=True):
    center = np.array(center)
    radius = np.linalg.norm(center - np.array(edge))

    total_sum, total_count, i, j = 0, 0, 0, 0
    rows, cols = img.shape

    while j < rows:
        while i < cols:
            if circle(i, j, center, radius) <= 0:
                total_sum += img[j, i]
                total_count += 1
            i += 1
        i = 0
        j += 1
    avg = total_sum / total_count

    voi_circle = patches.Circle(
        xy=center, radius=radius, linewidth=1, edgecolor="g", facecolor="none"
    )

    if display:
        fig, ax = plt.subplots()
        aplot = ax.imshow(img, cmap="gray")
        fig.colorbar(aplot, orientation="vertical")
        ax.add_patch(voi_circle)
        plt.show()

    return avg, voi_circle


def contrast_to_noise_tumor(
    image, voi_avg, center, edge, background_patch, display=True
):
    center = np.array(center)
    radius = np.linalg.norm(center - np.array(edge))
    print(center)
    print(radius)

    # get info for the background
    background_snr, background_avg, background_sigma = signal_to_noise(
        image, background_patch, display=False
    )

    print(f"Volume-of-Interest average: {voi_avg}")
    print(f"Background average: {background_avg}")
    print(f"Background standard deviation: {background_sigma}")
    cnr = (voi_avg - background_avg) / background_sigma
    print(f"contrast-to-noise ratio: {cnr}")

    if display:
        background_lowerleft, background_upperright = background_patch
        background_width, background_height = np.subtract(
            np.asarray(background_upperright), np.asarray(background_lowerleft)
        )

        fig, ax = plt.subplots()
        aplot = ax.imshow(image, cmap="gray")
        fig.colorbar(aplot, orientation="vertical")
        background_rect = patches.Rectangle(
            background_lowerleft,
            background_width,
            background_height,
            linewidth=1,
            edgecolor="r",
            facecolor="none",
        )
        ax.add_patch(background_rect)
        voi_circle = patches.Circle(
            xy=center, radius=radius, linewidth=1, edgecolor="g", facecolor="none"
        )
        ax.add_patch(voi_circle)
        plt.show()

    return cnr


source_file = "HGMTDerenzo"
vmax = None

################################################

# without tumor

# LAB 10,000th
# max pixel value = 78
# new max pixel value = 80
# source_file = r"/home/cameronpoe/Desktop/gauss_400px_1sig_z83.55_ta0.5_nola_tb_sphereatten.npy"
# vmax = 80

# LAB 1,000th
# max pixel value = 618
# new max pixel value = 630
# source_file = r"C:\Users\camer\Desktop\low-z_pet_scanner\brain_images\updated_3.14.23\without_tumor\e6_lab_gauss_400px_2sig_z83.55_nota_nola_sphereatten.npy"
# vmax = 630

# LAB 100th
# max pixel value = 5801
# new max pixel value = 6500
# source_file = r"C:\Users\camer\Desktop\low-z_pet_scanner\brain_images\updated_3.14.23\without_tumor\e7_lab_gauss_400px_2sig_z83.55_nota_nola_sphereatten.npy"
# vmax = 6500

################################################

# with tumor

# LAB 10,000th
# max pixel value = 115
# new max pixel value = 120
# source_file = r"/home/cameronpoe/Desktop/low-z_pet_scanner/brain_images/updated_3.14.23/with_tumor/e5_gauss_400px_2sig_z83.55_nota_nola_sphereatten.npy"
# vmax = 120

# LYSO 10,000th
# max pixel value = 24
# new max pixel value = 24
# source_file = r"/home/cameronpoe/Desktop/low-z_pet_scanner/brain_images/updated_3.14.23/updated_4.27.23_lyso_with_tumor/e5_gauss_400px_1sig_z83.55_ta0.5_nola_tb_sphereatten.npy"
# vmax = 24

# LAB 1,000th
# max pixel value = 916
# new max pixel value = 945
# source_file = r"C:\Users\camer\Desktop\low-z_pet_scanner\brain_images\updated_3.14.23\with_tumor\e6_gauss_400px_2sig_z83.55_nota_nola_sphereatten.npy"
# vmax = 945

# LYSO 1,000th
# max pixel value = 162
# new max pixel value = 172
# source_file = r"/home/cameronpoe/Desktop/low-z_pet_scanner/brain_images/updated_3.14.23/updated_4.27.23_lyso_with_tumor/e6_gauss_400px_1sig_z83.55_ta0.5_nola_tb_sphereatten.npy"
# vmax = 172

# LAB 100th
# max pixel value = 9126
# new max pixel value = 9750
source_file = r"output2.npy"
# vmax = 9750

# LYSO 100th
# max pixel value = 1382
# max pixel value = 1393
# source_file = r"/home/cameronpoe/Desktop/low-z_pet_scanner/brain_images/updated_3.14.23/updated_4.27.23_lyso_with_tumor/e7_gauss_400px_1sig_z83.55_ta0.5_nola_tb_sphereatten.npy"
# vmax = 1393

# 1/10,000th hgm brain, 1cm
# source_file = r'/home/cameronpoe/Desktop/hybrid_gamma_multiplier/hgm_lor_creator/brain_images/e5_1cm_100ps_gauss_400px_2sig_z83.55_nola_nota_xcatatten.npy'

# 1/1,0000th hgm brain, 1in
# source_file = r'/home/cameronpoe/Desktop/hybrid_gamma_multiplier/hgm_lor_creator/brain_images/e5_1in_100ps_gauss_400px_2sig_z83.55_nola_nota_xcatatten.npy'

# 1/1,000th hgm brain, 1cm
# source_file = r'/home/cameronpoe/Desktop/hybrid_gamma_multiplier/hgm_lor_creator/brain_images/e6_1cm_100ps_gauss_400px_2sig_z83.55_nola_nota_xcatatten.npy'

# 1/1,000th hgm brain, 1in
# source_file = r'/home/cameronpoe/Desktop/hybrid_gamma_multiplier/hgm_lor_creator/brain_images/e6_1in_100ps_gauss_400px_2sig_z83.55_nola_nota_xcatatten.npy'

# 1/100th hgm brain, 1cm
# source_file = r'/home/cameronpoe/Desktop/hybrid_gamma_multiplier/hgm_lor_creator/brain_images/e7_1cm_100ps_gauss_400px_2sig_z83.55_nola_nota_xcatatten.npy'

# 1/100th hgm brain, 1in
# source_file = r'/home/cameronpoe/Desktop/hybrid_gamma_multiplier/hgm_lor_creator/brain_images/e7_1in_100ps_gauss_400px_2sig_z83.55_nola_nota_xcatatten.npy'

############################################################################################################################

# 1/1,000th hgm brain, 1cm, 200ps FWHM
# source_file = r'/home/cameronpoe/Desktop/hybrid_gamma_multiplier/hgm_lor_creator/brain_images/200ps_images/1cm_gauss_400px_2sig_z83.55_nola_nota_xcatatten.npy'

# 1/1,000th hgm brain, 1in, 200ps FWHM
# source_file = r'/home/cameronpoe/Desktop/hybrid_gamma_multiplier/hgm_lor_creator/brain_images/200ps_images/1in_gauss_400px_2sig_z83.55_nola_nota_xcatatten.npy'

############################################################################################################################

# 1/1,000th hgm brain, 1cm, 10ns FWHM
# source_file = r'/home/cameronpoe/Desktop/hybrid_gamma_multiplier/hgm_lor_creator/brain_images/10ns_images/1cm_e6_gauss_400px_2sig_z83.55_nola_nota_xcatatten.npy'

# 1/1,000th hgm brain, 1in, 10ns FWHM
# source_file = r'/home/cameronpoe/Desktop/hybrid_gamma_multiplier/hgm_lor_creator/brain_images/10ns_images/1in_e6_gauss_400px_2sig_z83.55_nola_nota_xcatatten.npy'

################################################

# 1/1,000th hgm brain, 1in, 50ps FWHM
# source_file = r'/home/cameronpoe/Desktop/hybrid_gamma_multiplier/hgm_lor_creator/brain_images/50ps_images/1in_50ps_gauss_400px_2sig_z83.55_nola_nota_xcatatten.npy'

##############################################################################3

#### 100 ps
# 320um
# source_file = r'/home/cameronpoe/Desktop/hybrid_gamma_multiplier/hgm_lor_creator/brain_images/spatial_resolution_tests/100psFWHM_320umSigma_1inHGM_gauss_400px_2sig_z83.55_nola_nota_xcatatten.npy'

# 1 mm
# source_file = r'/home/cameronpoe/Desktop/hybrid_gamma_multiplier/hgm_lor_creator/brain_images/spatial_resolution_tests/100psFWHM_1mmSigma_1inHGM_gauss_400px_2sig_z83.55_nola_nota_xcatatten.npy'

# 3.2 mm
# source_file = r'/home/cameronpoe/Desktop/hybrid_gamma_multiplier/hgm_lor_creator/brain_images/spatial_resolution_tests/100psFWHM_3.2mmSigma_1inHGM_gauss_400px_2sig_z83.55_nola_nota_xcatatten.npy'

# 5 mm
# source_file = r'/home/cameronpoe/Desktop/hybrid_gamma_multiplier/hgm_lor_creator/brain_images/spatial_resolution_tests/100psFWHM_5mmSigma_1inHGM_gauss_400px_2sig_z83.55_nola_nota_xcatatten.npy'

#### 300 ps
# 320um
# source_file = r'/home/cameronpoe/Desktop/hybrid_gamma_multiplier/hgm_lor_creator/brain_images/spatial_resolution_tests/300psFWHM_320umSigma_1inHGM_gauss_400px_2sig_z83.55_nola_nota_xcatatten.npy'

# 1 mm
# source_file = r'/home/cameronpoe/Desktop/hybrid_gamma_multiplier/hgm_lor_creator/brain_images/spatial_resolution_tests/300psFWHM_1mmSigma_1inHGM_gauss_400px_2sig_z83.55_nola_nota_xcatatten.npy'

# 3.2 mm
# source_file = r'/home/cameronpoe/Desktop/hybrid_gamma_multiplier/hgm_lor_creator/brain_images/spatial_resolution_tests/300psFWHM_3.2mmSigma_1inHGM_gauss_400px_2sig_z83.55_nola_nota_xcatatten.npy'

# 5 mm
# source_file = r'/home/cameronpoe/Desktop/hybrid_gamma_multiplier/hgm_lor_creator/brain_images/spatial_resolution_tests/300psFWHM_5mmSigma_1inHGM_gauss_400px_2sig_z83.55_nola_nota_xcatatten.npy'

#### 500 ps
# 320 um
# source_file = r'/home/cameronpoe/Desktop/hybrid_gamma_multiplier/hgm_lor_creator/brain_images/spatial_resolution_tests/500psFWHM_320umSigma_1inHGM_gauss_400px_2sig_z83.55_nola_nota_xcatatten.npy'

# 1 mm
# source_file = r'/home/cameronpoe/Desktop/hybrid_gamma_multiplier/hgm_lor_creator/brain_images/spatial_resolution_tests/500psFWHM_1mmSigma_1inHGM_gauss_400px_2sig_z83.55_nola_nota_xcatatten.npy'

# 3.2 mm
# source_file = r'/home/cameronpoe/Desktop/hybrid_gamma_multiplier/hgm_lor_creator/brain_images/spatial_resolution_tests/500psFWHM_3.2mmSigma_1inHGM_gauss_400px_2sig_z83.55_nola_nota_xcatatten.npy'

# 5 mm
# source_file = r'/home/cameronpoe/Desktop/hybrid_gamma_multiplier/hgm_lor_creator/brain_images/spatial_resolution_tests/500psFWHM_5mmSigma_1inHGM_gauss_400px_2sig_z83.55_nola_nota_xcatatten.npy'


##############################
########## DERENZOS ##########
##############################

# 1/10,000th 1cm
# source_file = r'/home/cameronpoe/Desktop/hybrid_gamma_multiplier/hgm_lor_creator/derenzo_phantom/derenzo_images/e5_1cm_gauss_440px_2sig_z0_nota_nola_derenzoatten.npy'

# source_file = r'/home/cameronpoe/Desktop/hybrid_gamma_multiplier/hgm_lor_creator/derenzo_phantom/derenzo_images/e5_1in_gauss_440px_2sig_z0_nota_nola_derenzoatten.npy'

# 1/1,000th 1cm
# source_file = r'/home/cameronpoe/Desktop/hybrid_gamma_multiplier/hgm_lor_creator/derenzo_phantom/derenzo_images/e6_1cm_gauss_440px_2sig_z0_nota_nola_derenzoatten.npy'

# 1/1,000th 1in
# source_file = r'/home/cameronpoe/Desktop/hybrid_gamma_multiplier/hgm_lor_creator/derenzo_phantom/derenzo_images/e6_1in_gauss_440px_2sig_z0_nota_nola_derenzoatten.npy'

# 1/100th 1cm
# source_file = r'/home/cameronpoe/Desktop/hybrid_gamma_multiplier/hgm_lor_creator/derenzo_phantom/derenzo_images/e7_1cm_gauss_440px_2sig_z0_nota_nola_derenzoatten.npy'

# 1/100th 1in
# source_file = r'/home/cameronpoe/Desktop/hybrid_gamma_multiplier/hgm_lor_creator/derenzo_phantom/derenzo_images/e7_1in_gauss_440px_2sig_z0_nota_nola_derenzoatten.npy'

##############################


# source_file = r'/home/cameronpoe/Desktop/hybrid_gamma_multiplier/hgmt_sim_geant/LMCP_Full/data/pet_eff_tables/gauss_400px_2sig_z83.55_nola_nota_xcatatten.npy'
source_file = r"./HGMTDerenzo.npy"
pixels = 440
# filter_file = r"C:\Users\camer\topaspet-github\topas_truth_d\phantoms\xcat_topas_pet_examples\brain_images\filters\4_pi_filter.npy"
source_file_image = show_image(source_file, pixels, display=True, vmax=vmax)
# source_file_image_fft = fft(source_file_image, filter_file, display=False)
# source_file_image_fft_zerniked = plot_zernike_fitted(source_file_image, 3, 'z', 1, padding=75, display=False, padding_value=10)


# white_matter_patch = ([163, 185], [251, 151])
# gray_matter_patch = ([246, 224], [264, 207])

# avg, voi_circle = get_tumor_region(source_file_image, (184,273), (185, 254), display=False)
# contrast_to_noise_tumor(source_file_image, avg, (184,273), (185, 254), white_matter_patch, display=True)
# contrast_to_noise_tumor(source_file_image, avg, (184,273), (185, 254), gray_matter_patch, display=True)

# avg, voi_circle = get_tumor_region(source_file_image_fft_zerniked, (184,273), (185, 254), display=False)
# contrast_to_noise_tumor(source_file_image_fft_zerniked, avg, (184,273), (185, 254), white_matter_patch)


# white_snr, white_avg, white_sigma = signal_to_noise(source_file_image_fft_zerniked, white_matter_patch, display=True)
# contrast_to_noise(source_file_image_fft_zerniked, white_matter_patch, gray_matter_patch)
