# fits2comp
# Copyright (c) 2014-2022 Sebastiaan L. Zoutendijk
# Licensed under the terms of the MIT License, see file LICENSE

import argparse
from PIL import PngImagePlugin

from . import create_image, draw_annotations, mark_sources, read_colors, \
    transform_colors


def main():
    parser = argparse.ArgumentParser(
        description='Create composite color images from FITS files',
        epilog="""If only RED is given, GREEN and BLUE are taken equal to RED,
               which produces a grayscale image.  If RED and GREEN are given,
               BLUE is extrapolated from them.""")
    parser.add_argument('-A', '--no-annotation-area', action='store_true',
                        help='draw ruler and compass in the image area')
    parser.add_argument('-b', '--background', default=0, type=float,
                        help='background pixel value in FITS')
    parser.add_argument('-C', '--compass', action='store_true',
                        help='add compass indicating north and east')
    parser.add_argument('-c', '--cutoff', default=65535, type=float,
                        help='cut-off pixel value in FITS')
    parser.add_argument('-l', '--source-list', metavar='FILENAME', type=str,
                        help='indicate source positions from file')
    parser.add_argument('-n', '--name', metavar='NAME', type=str,
                        help='add target name')
    parser.add_argument('-o', '--offset', metavar=('X', 'Y'), nargs=2,
                        default=(0, 0), type=float,
                        help='X and Y offset of source positions, in pixels')
    parser.add_argument('-R', '--ruler', metavar=('ANGULAR', 'PHYSICAL'),
                        nargs=2, type=str,
                        help='add ruler of ANGULAR length, with ANGULAR ' + \
                            'and PHYSICAL labels')
    parser.add_argument('-S', '--stretch', default=0, type=float,
                        help='stretch around background with this amplitude')
    parser.add_argument('-s', '--softening', default=1, type=float,
                        help='softening parameter')
    parser.add_argument('red', metavar='RED', help='red channel FITS file')
    parser.add_argument('green', metavar='GREEN', nargs='?',
                        help='green channel FITS file')
    parser.add_argument('blue', metavar='BLUE', nargs='?',
                        help='blue channel FITS file')
    parser.add_argument('composite', metavar='COMPOSITE',
                        help='output composite PNG image')
    sys_args = parser.parse_args()

    # Save the arguments as PNG tags for reference
    tags = PngImagePlugin.PngInfo()
    tags.add_text('Software', 'fits2comp')
    tags.add_text('--background', repr(sys_args.background))
    tags.add_text('--compass', repr(sys_args.compass))
    tags.add_text('--cutoff', repr(sys_args.cutoff))
    tags.add_text('--name', repr(sys_args.name))
    tags.add_text('--no-annotation-area', repr(sys_args.no_annotation_area))
    tags.add_text('--offset', repr(sys_args.offset))
    tags.add_text('--ruler', repr(sys_args.ruler))
    tags.add_text('--stretch', repr(sys_args.stretch))
    tags.add_text('--softening', repr(sys_args.softening))
    tags.add_text('--source-list', repr(sys_args.source_list))
    tags.add_text('RED', repr(sys_args.red))
    tags.add_text('BLUE', repr(sys_args.blue))
    tags.add_text('GREEN', repr(sys_args.green))
    tags.add_text('COMPOSITE', repr(sys_args.composite))

    red, green, blue, mask, wcs = \
        read_colors(sys_args.red, sys_args.green, sys_args.blue,
                    sys_args.background)
    Red, Green, Blue = \
        transform_colors(red, green, blue, sys_args.cutoff, sys_args.stretch,
                         sys_args.softening)
    image = create_image(Red, Green, Blue, mask)
    image = \
        mark_sources(image, mask, wcs, sys_args.source_list, sys_args.offset)
    image = \
        draw_annotations(image, wcs, sys_args.name, sys_args.ruler,
                         sys_args.compass, sys_args.no_annotation_area)

    image.save(sys_args.composite, format='PNG', pnginfo=tags)


if __name__ == '__main__':
    main()
