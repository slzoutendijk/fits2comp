# fits2comp
# Copyright (c) 2014-2022 Sebastiaan L. Zoutendijk
# Licensed under the terms of the MIT License, see file LICENSE

from astropy.coordinates import Angle, SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
import numpy as np
from PIL import Image, ImageDraw, ImageFont


def create_image(Red, Green, Blue, Alpha):
    # Take red, green, blue, and transparency channels [0, 1] and create an
    # RGBA image

    Red_image = Image.fromarray(np.around(255*Red[::-1,:]).astype('uint8'))
    Green_image = \
        Image.fromarray(np.around(255*Green[::-1,:]).astype('uint8'))
    Blue_image = Image.fromarray(np.around(255*Blue[::-1,:]).astype('uint8'))
    Alpha_image = \
        Image.fromarray(np.around(255*Alpha[::-1,:]).astype('uint8'))
    image = \
        Image.merge('RGBA', (Red_image, Green_image, Blue_image, Alpha_image))

    return image


def draw_annotations(image, wcs, name, ruler, compass, no_annotation_area):
    # Create an annotation area and draw a ruler, compass, and object name, if
    # necessary
    if ((name is not None or ruler is not None or compass)
            and not no_annotation_area):
        image = make_annotation_area(image)

    if ruler is not None:
        image, left = draw_ruler(image, wcs, ruler)
    else:
        left = 0 # start of object name area

    if compass:
        image, right = draw_compass(image, wcs)
    else:
        right = image.size[0] # end of object name area

    if name is not None:
        image = draw_name(image, name, loc=np.mean([left, right]))

    return image


def draw_compass(image, wcs):
    # Draw a compass indicating north and east

    draw = ImageDraw.Draw(image)

    # Center of compass
    compass_start_px = (image.size[0]-30, 30)
    compass_start_deg = \
        SkyCoord(*wcs.all_pix2world(*compass_start_px, 0), unit='deg')
    compass_left_px = (image.size[0]-60, 30)

    # North arm of compass
    compass_north_test_deg = \
        compass_start_deg.directional_offset_by(0, Angle('1arcsec'))
    compass_north_test_px = \
        wcs.all_world2pix(compass_north_test_deg.ra.deg,
                          compass_north_test_deg.dec.deg, 0)
    compass_north_xy_ang = \
        np.arctan2(compass_north_test_px[1]-compass_start_px[1],
                   compass_north_test_px[0]-compass_start_px[0])
    compass_north_end_px = \
        (compass_start_px[0] + 20 * np.cos(compass_north_xy_ang),
         compass_start_px[1] + 20 * np.sin(compass_north_xy_ang))
    compass_north_text_px = \
        (compass_start_px[0] + \
             15*np.sqrt(2) * np.cos(compass_north_xy_ang-np.pi/4),
         compass_start_px[1] + \
             15*np.sqrt(2) * np.sin(compass_north_xy_ang-np.pi/4))

    # East arm of compass
    compass_east_test_deg = \
        compass_start_deg.directional_offset_by(np.pi/2, Angle('1arcsec'))
    compass_east_test_px = \
        wcs.all_world2pix(compass_east_test_deg.ra.deg,
                          compass_east_test_deg.dec.deg, 0)
    compass_east_xy_ang = \
        np.arctan2(compass_east_test_px[1]-compass_start_px[1],
                   compass_east_test_px[0]-compass_start_px[0])
    compass_east_end_px = \
        (compass_start_px[0] + 20 * np.cos(compass_east_xy_ang),
         compass_start_px[1] + 20 * np.sin(compass_east_xy_ang))
    compass_east_text_px = \
        (compass_start_px[0] + \
            15*np.sqrt(2) * np.cos(compass_east_xy_ang+np.pi/4),
         compass_start_px[1] + \
            15*np.sqrt(2) * np.sin(compass_east_xy_ang+np.pi/4))

    # Draw arms, arrow heads, and labels
    draw.line(((compass_north_end_px[0],
                image.size[1]-compass_north_end_px[1]),
               (compass_start_px[0], image.size[1]-compass_start_px[1]),
               (compass_east_end_px[0], image.size[1]-compass_east_end_px[1])),
              fill='black', width=5)
    draw.regular_polygon(((compass_north_end_px[0],
                           image.size[1]-compass_north_end_px[1]),
                          5),
                         3,
                         compass_north_xy_ang*180/np.pi-90,
                         fill='black')
    draw.regular_polygon(((compass_east_end_px[0],
                           image.size[1]-compass_east_end_px[1]),
                          5),
                         3,
                         compass_east_xy_ang*180/np.pi-90,
                         fill='black')
    draw.text((compass_north_text_px[0],
               image.size[1]-compass_north_text_px[1]),
              'N',
              font=ImageFont.truetype('DejaVuSans.ttf', 24),
              fill='black',
              anchor='mm')
    draw.text((compass_east_text_px[0], image.size[1]-compass_east_text_px[1]),
              'E', font=ImageFont.truetype('DejaVuSans.ttf', 24),
              fill='black', anchor='mm')

    return image, compass_left_px[0]

def draw_name(image, name, loc=None):
    # Draw the object name

    draw = ImageDraw.Draw(image)
    if loc is None:
        name_px = (image.size[0]/2, 30)
    else:
        name_px = (loc, 30)
    draw.text((name_px[0], image.size[1]-name_px[1]), name,
              font=ImageFont.truetype('DejaVuSans.ttf', 24), fill='black',
              anchor='mm')

    return image


def draw_ruler(image, wcs, ruler):
    # Draw a ruler indicating angular and physical scale

    draw = ImageDraw.Draw(image)

    line_start_px = (0, 30)
    line_start_deg = \
        SkyCoord(*wcs.all_pix2world(*line_start_px, 0), unit='deg')
    line_test_px = (1, 30)
    line_test_deg = SkyCoord(*wcs.all_pix2world(*line_test_px, 0), unit='deg')
    line_pos_ang = line_start_deg.position_angle(line_test_deg)
    line_end_deg = \
        line_start_deg.directional_offset_by(line_pos_ang, Angle(ruler[0]))
    line_end_px = \
        wcs.all_world2pix(line_end_deg.ra.deg, line_end_deg.dec.deg, 0)
    draw.line(((line_start_px[0], image.size[1]-line_start_px[1]),
               (line_end_px[0], image.size[1]-line_end_px[1])),
              fill='black',
              width=5)

    half_line = np.mean([line_start_px, line_end_px], axis=0)
    draw.text((half_line[0], image.size[1]-half_line[1]-15), ruler[0],
              font=ImageFont.truetype('DejaVuSans.ttf', 24), fill='black',
              anchor='mm')
    draw.text((half_line[0], image.size[1]-half_line[1]+15), ruler[1],
              font=ImageFont.truetype('DejaVuSans.ttf', 24), fill='black',
              anchor='mm')

    return image, line_end_px[0]


def map_colors(red, green, blue, cutoff, softening):
    # Map the colors to [0, 1] using the Lupton et al. algorithm

    Red = np.zeros(red.shape)
    Green = np.zeros(green.shape)
    Blue = np.zeros(blue.shape)

    average = np.mean([red, green, blue], axis=0)
    pos = (average > 0)

    mapped = mapping(average[pos], cutoff, softening)
    Red[pos] = red[pos] * mapped / average[pos]
    Green[pos] = green[pos] * mapped / average[pos]
    Blue[pos] = blue[pos] * mapped / average[pos]

    # Clip to [0, 1], but preserve the color at the bright end
    amax = np.amax([Red, Green, Blue], axis=0)
    Red[amax > 1] /= amax[amax > 1]
    Green[amax > 1] /= amax[amax > 1]
    Blue[amax > 1] /= amax[amax > 1]
    Red[Red < 0] = 0
    Green[Green < 0] = 0
    Blue[Blue < 0] = 0

    return Red, Green, Blue


def make_annotation_area(image):
    # Expand the image at the bottom to accommodate annotations

    image = image.crop((0, 0, image.size[0], image.size[1]+60))

    return image


def mapping(value, cutoff, softening):
    # Map the intensity following Lupton et al.

    mapped = \
        scaling(value.astype(float), softening) / scaling(cutoff, softening)

    mapped[mapped < 0] = 0
    mapped[mapped > 1] = 1

    return mapped


def mark_sources(image, mask, wcs, source_list, offset):
    # Mark sources from a list with circles

    if source_list is not None:
        draw = ImageDraw.Draw(image)
        cat = Table.read(source_list) # anything tabular goes
        sources = \
            wcs.all_world2pix(
                np.array([cat['RAdeg'], cat['DEdeg']]).transpose(),
                0)

        for source in sources:
            source = np.array([source[0]+offset[0], (source[1]+offset[1])])
            for x in range(mask.shape[1]):
                for y in range(mask.shape[0]):
                    # Mark only if there is unmasked data inside the circle
                    if ((x-source[0])**2 + (y-source[1])**2 <= 10**2
                            and mask[y,x]):
                        source = np.array([source[0], image.size[1]-source[1]])
                        draw.ellipse((tuple(source-10), tuple(source+10)),
                                     outline='yellow', width=1)
                        break # data found, jump to (A)
                else: # no data found yet at this x for any y
                    continue # advance x
                # (A)
                break # go to next source

    return image


def read_color(fname, background, mask=None, return_wcs=False):
    # Read a FITS file, assuming the relevant extension is called DATA

    color_file = fits.open(fname)

    ext = color_file.index_of('DATA')
    if return_wcs:
        wcs = WCS(color_file[ext].header)
    if mask is None:
        mask = np.isfinite(color_file[ext].data)
    else:
        mask &= np.isfinite(color_file[ext].data)
    color = np.nan_to_num(color_file[ext].data) - background

    color_file.close()

    if return_wcs:
        return color, mask, wcs
    else:
        return color, mask


def read_colors(red_name, green_name, blue_name, background):
    # Read all the colors, and extrapolate for those not given

    red, mask, wcs = read_color(red_name, background, return_wcs=True)
    if green_name is not None:
        green, mask = read_color(green_name, background, mask=mask)
        if blue_name is not None:
            blue, mask = read_color(blue_name, background, mask=mask)
        else:
            blue = 2 * green - red
    else:
        green = red
        blue = red

    return red, green, blue, mask, wcs


def scaling(value, softening):
    # Scaling function of Lupton et al.

    scaled = np.arcsinh(value/softening)

    return scaled


def transform_colors(red, green, blue, cutoff, stretch, softening):
    # Given colors, transform these to [0, 1] using the Lupton et al. mappings

    if stretch == 0:
        Red, Green, Blue = map_colors(red, green, blue, cutoff, softening)
    else:
        # Place the background level at 0.5 and symmetrically stretch around
        # this level
        posRed, posGreen, posBlue = \
            map_colors(+red, +green, +blue, stretch, softening)
        negRed, negGreen, negBlue = \
            map_colors(-red, -green, -blue, stretch, softening)
        Red, Green, Blue = \
            0.5 + \
            np.mean([[posRed, posGreen, posBlue], [negRed, negGreen, negBlue]],
                    axis=0)

    return Red, Green, Blue


