# fits2comp

fits2comp creates composite color images from FITS files using the algorithm of
Lupton et al. ([2004][L04]).  Though the algorithm has been implemented in
Astropy since I wrote the first version of this program,  the main benefit of
fits2comp is that you can call it with a number of options that will take care
of a number of things for you -- no need to code things up.

[L04]: https://ui.adsabs.harvard.edu/abs/2004PASP..116..133L/abstract

You will need to create FITS images in up to three filters by yourself.  The
FITS images need to be on the same grid and are combined pixel by pixel, WCS
information is ignored for this purpose.  The WCS of the red channel is used
for more advanced functionality: drawing a ruler and a compass, and marking
source positions from a file readable as an Astropy table.  For my own science
with the MUSE integral-field spectrograph, I use the `exp_combine` recipe or
the command `muse_cube_filter` to create images in several filters from a data
cube.

fits2comp is a Python 3 package, installable via `pip install fits2comp`.  It
requires Astropy, NumPy, and Pillow-PIL, which will be automatically installed.
It uses argparse, so it accepts the `--help` option.  The most important
options are `-c|--cutoff` and `-s|--softening`, which together set the color
scale of the output image.  You will have to try a few values to get an image
that you like.

Pre-released versions of this package were used to create the images in
Zoutendijk et al. ([2020][PaperI], Fig. 1; [2021a][PaperII], Fig. 1;
[2021b][PaperIII], Fig. 1) As an example, I created the Antlia B image in
Zoutendijk et al. ([2021b][PaperIII], Fig. 1) with

```
fits2comp \
    --cutoff 10 --softening 1 \
    --ruler "30 arcsec" "196 pc" --name "Ant B" --compass \
    --source-list AntB.csv --offset 5 5 \
    AntB_im/DATACUBE_exp_combine-nomask-noautocal_SDSS_i.fits \
    AntB_im/DATACUBE_exp_combine-nomask-noautocal_SDSS_r.fits \
    AntB_im/DATACUBE_exp_combine-nomask-noautocal_SDSS_g.fits \
    AntB_gri.png
```

[PaperI]: https://ui.adsabs.harvard.edu/abs/2020A%26A...635A.107Z/abstract
[PaperII]: https://ui.adsabs.harvard.edu/abs/2021A%26A...651A..80Z/abstract
[PaperIII]: https://ui.adsabs.harvard.edu/abs/2021arXiv211209374Z/abstract

You can find the options used to create a composite image in the
image metadata, in case you forget what options you used after
painstakingly finding the optimal combination for your target.

## Contact

In case of bugs or questions, please file an issue on GitHub or
send me an e-mail: `zoutendijk@strw.leidenuniv.nl`.

## License

This project is licensed under the terms of the [MIT License](LICENSE).
