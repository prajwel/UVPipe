#!/usr/bin/env python3

import re
import os
import fnmatch

from glob import glob
from astropy.io import fits
from shutil import copy, rmtree, move


def finD(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in dirs:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result


def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result


# Place the driver module param file.
L1_filename = glob("LEVL1AS1UVT*tar_V*.2*")[0]
L1_datepid = list(L1_filename[:26])
L1_datepid[4] = "2"
ICD_driver_name = "".join(L1_datepid) + "_L2_DM_params.txt"
os.rename("UVIT_DriverModule.par", ICD_driver_name)
os.mkdir("pipeline")
move(ICD_driver_name, "./pipeline/")

# To find output FUV & NUV directories.
fuvdlist = glob("output_FUV_[0-9]*")
nuvdlist = glob("output_NUV_[0-9]*")
if len(fuvdlist) == 0 and len(fuvdlist) == 0:
    print('\nCheck if output files in the format "output_[F,N]UV*" are present\n')

# To make a dictionary of all the VIS RAS files.
rasdlist = find("*dr.fits", ".")
ras_dict = dict([(rasd[-75:-68], rasd[1:]) for rasd in rasdlist])

# To find and rename the orbit-wise FUV, NUV files
currdir = os.getcwd()
for fuvd in fuvdlist:
    os.chdir(currdir + "/" + fuvd)
    Fasi = find("*as_Sig.fits", ".")
    Fsnri = find("*l2_radec.fits", ".")
    Fsigi = find("*sigFlipped_rg.fits", ".")
    FAexpi = find("*expFlipped_rg.fits", ".")
    FAerri = find("*as_NoiseMap.fits", ".")

    if len(Fsnri) == 1:
        move(Fsnri[0], ".")

    if len(FAexpi) == 1:
        move(FAexpi[0], ".")

    if len(Fasi) == 1:
        move(Fasi[0], ".")

    if len(FAerri) == 1:
        move(FAerri[0], ".")

    if len(Fsigi) == 1:
        move(Fsigi[0], ".")
        fuvdir = glob("uvit")
        if len(fuvdir) != 0:
            rmtree(fuvdir[0])

    if len(Fasi) == 0 and len(Fsigi) > 0:
        print("No frame with astrometry in {}".format(fuvd))

    if len(Fsigi) == 0:
        fuvdir = glob("uvit")
        if len(fuvdir) != 0:
            print("directory {} is empty".format(fuvd))
            rmtree(fuvdir[0])
            continue

    if len(Fsigi) > 1:
        print("\nExists more than one FUV image inside {}".format(fuvd))
        print("Formatting not carried out.\n")
        continue

    Fas = glob("*as_Sig.fits")
    Fsnr = glob("*l2_radec.fits")
    Fsig = glob("*sigFlipped_rg.fits")
    FAexp = glob("*expFlipped_rg.fits")
    FAerr = glob("*as_NoiseMap.fits")

    # To get the corresponding RAS file.
    if len(Fsnr) == 1:
        Fhdu = fits.open(Fsnr[0])
        history = ""
        if "history" in Fhdu[0].header:
            history = Fhdu[0].header["history"]
        else:
            print("History missing for events list in {}".format(fuvd))
            if len(Fas) == 1:
                Fhdu = fits.open(Fas[0])
                history = Fhdu[0].header["history"]
            else:
                print("No way to get RAS in {}\n".format(fuvd))

        if len(history) > 0:
            re_result = re.search(r"V/uvtV\.\d{2}/uvtC", str(history))

            if type(re_result.group()) is str:
                ras_num = re_result.group()[2:9]
            else:
                print("\nmultiple RAS?! Big problem!")
                exit()

            ras_file = currdir + ras_dict[ras_num]
            copy(ras_file, ".")
    else:
        print("\nFUV events file is absent at {}!\n".format(fuvd))
        continue

    ras = glob("*dr.fits")

    if len(Fas) != 0:
        os.rename(Fas[0], (Fas[0][0:36] + "A_l2img.fits"))
    if len(Fsnr) != 0:
        os.rename(Fsnr[0], (Fsnr[0][0:36] + "_l2ce.fits"))
    if len(Fsig) != 0:
        os.rename(Fsig[0], (Fsig[0][0:36] + "I_l2img.fits"))
    if len(FAexp) != 0:
        os.rename(FAexp[0], (FAexp[0][0:36] + "I_l2exp.fits"))
    if len(FAerr) != 0:
        os.rename(FAerr[0], (FAerr[0][0:36] + "A_l2err.fits"))
    if len(ras) != 0:
        os.rename(ras[0], (ras[0][0:36] + "_l2dr.fits"))


os.chdir(currdir)
for nuvd in nuvdlist:
    os.chdir(currdir + "/" + nuvd)
    Nasi = find("*as_Sig.fits", ".")
    Nsnri = find("*l2_radec.fits", ".")
    Nsigi = find("*sigFlipped_rg.fits", ".")
    NAexpi = find("*expFlipped_rg.fits", ".")
    NAerri = find("*as_NoiseMap.fits", ".")

    if len(Nsnri) == 1:
        move(Nsnri[0], ".")

    if len(NAexpi) == 1:
        move(NAexpi[0], ".")

    if len(Nasi) == 1:
        move(Nasi[0], ".")

    if len(NAerri) == 1:
        move(NAerri[0], ".")

    if len(Nsigi) == 1:
        move(Nsigi[0], ".")
        nuvdir = glob("uvit")
        if len(nuvdir) != 0:
            rmtree(nuvdir[0])

    if len(Nasi) == 0 and len(Nsigi) > 0:
        print("No frame with astrometry in {}".format(nuvd))

    if len(Nsigi) == 0:
        nuvdir = glob("uvit")
        if len(nuvdir) != 0:
            print("directory {} is empty".format(nuvd))
            rmtree(nuvdir[0])
            continue

    if len(Nsigi) > 1:
        print("\nExists more than one NUV image inside {}".format(nuvd))
        print("Formatting not carried out.\n")
        continue

    Nas = glob("*as_Sig.fits")
    Nsnr = glob("*l2_radec.fits")
    Nsig = glob("*sigFlipped_rg.fits")
    NAexp = glob("*expFlipped_rg.fits")
    NAerr = glob("*as_NoiseMap.fits")

    # To get the corresponding RAS file.
    if len(Nsnr) == 1:
        Nhdu = fits.open(Nsnr[0])
        history = ""
        if "history" in Nhdu[0].header:
            history = Nhdu[0].header["history"]
        else:
            print("History missing for events list in {}".format(nuvd))
            if len(Nas) == 1:
                Nhdu = fits.open(Nas[0])
                history = Nhdu[0].header["history"]
            else:
                print("No way to get RAS in {}\n".format(nuvd))

        if len(history) > 0:
            re_result = re.search(r"V/uvtV\.\d{2}/uvtC", str(history))

            if type(re_result.group()) is str:
                ras_num = re_result.group()[2:9]
            else:
                print("\nmultiple RAS?! Big problem!")
                exit()

            ras_file = currdir + ras_dict[ras_num]
            copy(ras_file, ".")
    else:
        print("\nNUV events file is absent at {}!\n".format(nuvd))
        continue

    ras = glob("*dr.fits")

    if len(Nas) != 0:
        os.rename(Nas[0], (Nas[0][0:36] + "A_l2img.fits"))
    if len(Nsnr) != 0:
        os.rename(Nsnr[0], (Nsnr[0][0:36] + "_l2ce.fits"))
    if len(Nsig) != 0:
        os.rename(Nsig[0], (Nsig[0][0:36] + "I_l2img.fits"))
    if len(NAexp) != 0:
        os.rename(NAexp[0], (NAexp[0][0:36] + "I_l2exp.fits"))
    if len(NAerr) != 0:
        os.rename(NAerr[0], (NAerr[0][0:36] + "A_l2err.fits"))
    if len(ras) != 0:
        os.rename(ras[0], (ras[0][0:36] + "_l2dr.fits"))

os.chdir(currdir)

# renaming the FUV & NUV folders.
for fuvd in fuvdlist:
    try:
        t = int(fuvd[-2:])
        os.rename(fuvd, "F_" + str(t))
    except ValueError:
        os.rename(fuvd, "F_0" + fuvd[-1:])

for nuvd in nuvdlist:
    try:
        t = int(nuvd[-2:])
        os.rename(nuvd, "N_" + str(t))
    except ValueError:
        os.rename(nuvd, "N_0" + nuvd[-1:])

# Creating structure as per ICD
F_dirs = glob("F_*[0-9]")
N_dirs = glob("N_*[0-9]")
F_dir_max_index = max([F_dir[-2:] for F_dir in F_dirs])
N_dir_max_index = max([N_dir[-2:] for N_dir in N_dirs])
fcount = max(F_dir_max_index, N_dir_max_index)
fcount = int(fcount)

for fc in range(1, fcount + 1):
    if fc < 10:
        try:
            uvtD = "uvt_0" + str(fc)
            VD = "V_0" + str(fc)
            os.mkdir(uvtD)
            os.mkdir(VD)
        except OSError:
            print("Check user permissions")

        try:
            move("F_0" + str(fc), uvtD)
        except IOError:
            pass

        try:
            move("N_0" + str(fc), uvtD)
        except IOError:
            pass

        for dras in find("*dr.fits", "./" + uvtD):
            copy(dras, VD)
            os.remove(dras)

        move(VD, uvtD)

    else:
        try:
            uvtD = "uvt_" + str(fc)
            VD = "V_" + str(fc)
            os.mkdir(uvtD)
            os.mkdir(VD)
        except OSError:
            print("Check user permissions")

        try:
            move("F_" + str(fc), uvtD)
        except IOError:
            pass

        try:
            move("N_" + str(fc), uvtD)
        except IOError:
            pass

        for dras in find("*dr.fits", "./" + uvtD):
            copy(dras, VD)
            os.remove(dras)

        move(VD, uvtD)

# To find, move, and rename the combined FUV, NUV files.
os.mkdir("uvt_ci/")
os.mkdir("uvt_ci/F_ci/")
os.mkdir("uvt_ci/N_ci/")
for ciD in finD("*outputFUV*", "."):
    os.chdir(currdir + ciD[1:])
    cis = find("*.fits", ".")
    for ci in cis:
        move(ci, currdir + "/uvt_ci/F_ci/")

os.chdir(currdir)

for ciD in finD("*outputNUV*", "."):
    os.chdir(currdir + ciD[1:])
    cis = find("*.fits", ".")
    for ci in cis:
        move(ci, currdir + "/uvt_ci/N_ci/")

os.chdir(currdir + "/uvt_ci/")

for ci_as in find("*as_Sig.fits", "."):
    as_hdu = fits.open(ci_as)
    prefx = as_hdu[0].header["NAMEPRFX"]
    os.rename(ci_as, ci_as[:7] + prefx[:36] + "A_l2imb.fits")

for ci_exp in find("*FinalImage_Exp.fits", "."):
    os.remove(ci_exp)

for ci_err in find("*FinalImage_NoiseMap.fits", "."):
    os.remove(ci_err)

for ci_sig in find("*FinalImage_Sig.fits", "."):
    sig_hdu = fits.open(ci_sig)
    prefx = sig_hdu[0].header["NAMEPRFX"]
    os.rename(ci_sig, ci_sig[:7] + prefx[:36] + "A_l2ima.fits")

for ci_as_exp in find("*as_Exp.fits", "."):
    exp_hdu = fits.open(ci_as_exp)
    prefx = exp_hdu[0].header["NAMEPRFX"]
    os.rename(ci_as_exp, ci_as_exp[:7] + prefx[:36] + "A_l2exb.fits")

for ci_as_err in find("*as_NoiseMap.fits", "."):
    err_hdu = fits.open(ci_as_err)
    prefx = err_hdu[0].header["NAMEPRFX"]
    os.rename(ci_as_err, ci_as_err[:7] + prefx[:36] + "A_l2erb.fits")


os.chdir(currdir)
rm1dir = glob("outputFUV*")
rm2dir = glob("outputNUV*")
for rm1 in rm1dir:
    rmtree(rm1)
for rm2 in rm2dir:
    rmtree(rm2)

for uvt_dir in glob("uvt_[0-9]*"):
    dir_num = uvt_dir.split("uvt_")[-1]
    fuv_dir = f"{uvt_dir}/F_{dir_num}"
    if not os.path.exists(fuv_dir):
        os.makedirs(fuv_dir)
    nuv_dir = f"{uvt_dir}/N_{dir_num}"
    if not os.path.exists(nuv_dir):
        os.makedirs(nuv_dir)
    vis_dir = f"{uvt_dir}/V_{dir_num}"
    if not os.path.exists(vis_dir):
        os.makedirs(vis_dir)

print("\nDone.")
