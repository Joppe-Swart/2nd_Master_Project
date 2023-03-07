"""

Author: Joppe Swart
Created: January 2023
Last modified: January 2023
Description: This scrip is meant to renama all files so we have one long 

"""

import os


def rename(sourcename, channels):
    """
    rename: Function that renames the filenames of split data blocks
            The funtion numbers two files at the same time.
    INPUTS:
        sourcename: the name of the input block
        channels: Number of channels the input dir contains.
    """


    print(f'Renaming the files')
    if sourcename=='BC_I':
        teller = 0
        for i in range(channels):
            if teller == channels:
                break
            print(f'Rename channel {teller} and {teller+channels}')
            if not os.path.exists(f"{teller:04d}-Q-image.fits"):
                os.rename(f"{sourcename}-Split-1-{teller:04d}-Q-image.fits",\
                          f"{teller:04d}-Q-image.fits")
            if not os.path.exists(f"{teller+channels:04d}-Q-image.fits"):
                os.rename(f"{sourcename}-Split-2-{teller:04d}-Q-image.fits",\
                          f"{teller+channels:04d}-Q-image.fits")
            if not os.path.exists(f"{teller:04d}-U-image.fits"):
                os.rename(f"{sourcename}-Split-1-{teller:04d}-U-image.fits",\
                          f"{teller:04d}-U-image.fits")
            if not os.path.exists(f"{teller+channels:04d}-U-image.fits"):
                os.rename(f"{sourcename}-Split-2-{teller:04d}-U-image.fits",\
                          f"{teller+channels:04d}-U-image.fits")
            if not os.path.exists(f"{teller:04d}-I-image.fits"):
                os.rename(f"{sourcename}-Split-1-{teller:04d}-image.fits",\
                          f"{teller:04d}-I-image.fits")
            if not os.path.exists(f"{teller+channels:04d}-I-image.fits"):
                os.rename(f"{sourcename}-Split-2-{teller:04d}-image.fits",\
                          f"{teller+channels:04d}-I-image.fits")
            teller +=1

    if sourcename=='BC_II':
        teller = 0
        teller_all=103
        for i in range(channels):
            if teller == channels:
                break
            print(f'Rename channel {teller} and {teller+channels}')
            teller_all +=1
            if not os.path.exists(f"{teller_all:04d}-Q-image.fits"):
                os.rename(f"{sourcename}_Split_1-{teller:04d}-Q-image.fits",\
                          f"{teller_all:04d}-Q-image.fits")
            if not os.path.exists(f"{teller_all+channels:04d}-Q-image.fits"):
                os.rename(f"{sourcename}_Split_2-{teller:04d}-Q-image.fits",\
                          f"{teller_all+channels:04d}-Q-image.fits")
            if not os.path.exists(f"{teller_all:04d}-U-image.fits"):
                os.rename(f"{sourcename}_Split_1-{teller:04d}-U-image.fits",\
                          f"{teller_all:04d}-U-image.fits")
            if not os.path.exists(f"{teller_all+channels:04d}-U-image.fits"):
                os.rename(f"{sourcename}_Split_2-{teller:04d}-U-image.fits",\
                          f"{teller_all+channels:04d}-U-image.fits")
            if not os.path.exists(f"{teller_all:04d}-I-image.fits"):
                os.rename(f"{sourcename}_Split_1-{teller:04d}-image.fits",\
                          f"{teller_all:04d}-I-image.fits")
            if not os.path.exists(f"{teller_all+channels:04d}-I-image.fits"):
                os.rename(f"{sourcename}_Split_2-{teller:04d}-image.fits",\
                          f"{teller_all+channels:04d}-I-image.fits")
            teller+=1

    if sourcename=='BC_III':
        teller = 0
        teller_all=227
        for i in range(channels):
            if teller == channels:
                break
            print(f'Rename channel {teller} and {teller+channels}')
            teller_all +=1
            if not os.path.exists(f"{teller_all:04d}-Q-image.fits"):
                os.rename(f"{sourcename}-{teller:04d}-Q-image.fits",\
                          f"{teller_all:04d}-Q-image.fits")
            if not os.path.exists(f"{teller_all:04d}-U-image.fits"):
                os.rename(f"{sourcename}-{teller:04d}-U-image.fits",\
                          f"{teller_all:04d}-U-image.fits")
            if not os.path.exists(f"{teller_all:04d}-I-image.fits"):
                os.rename(f"{sourcename}-{teller:04d}-image.fits",\
                          f"{teller_all:04d}-I-image.fits")
            teller+=1

rename(sourcename='BC_I', channels=52)
rename(sourcename='BC_II', channels=62)
rename(sourcename='BC_III', channels=25)

