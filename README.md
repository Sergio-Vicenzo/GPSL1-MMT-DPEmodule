# MMT-integrated DPE
An open-source GNSS Multipath Mitigation Technology-integrated Direct Position Estimation Plug-in Module for Two-Step Positioning SDRs, integrated into SoftGNSS v3.0 by Borre et al. (2007)

[Intelligent Positioning and Navigation Laboratory (IPNL)](https://www.polyu.edu.hk/aae/ipn-lab/us/index.html) / [PNT Signal Processing Group](https://pbingxu.github.io/team/)

The Hong Kong Polytechnic University

## Related Publication
---

## Introduction
MMT-integrated DPE is an extension of our previously introduced [DPE_module v1.0](https://github.com/Sergio-Vicenzo/GPSL1-DPEmodule).

Though DPE has been proven to be robust against MP, previous research has proved that its superior performance against two-step positioning (2SP) typically falters in deep urban environments where multipath (MP) and non-line-of-sight (NLOS) are the majority of signals received (Vicenzo et al. 2023). This does not necessarily mean that the performance of DPE is worse than 2SP, but rather that its performance is depreciated to a large degree with increasing errors from MP and NLOS the same way 2SP does. 

Tang et al. (2024) had also recently showed that while DPE remains more accurate in comparison to 2SP in harsh cases such as 4 out of 8 satellites being MP, DPE error still reaches up to tens, or even hundreds of meters with NLOS measurements, which makes it definitely unsuitable for urban positioning. With research into mitigating MP with 2SP has been widely explored, the prevalence and need for DPE as a robust positioning method against MP has diminished significantly, especially with its high computational load and inapplicability to commercial receivers producing RINEX-level measurements.

To solve this issue for DPE, we to introduce a Multipath Mitigation Technology (MMT)-integrated DPE. MMT was introduced as an efficient estimator for accurate estimation of the code delays and carrier phase of the LOS and reflected signal (Weill 2002). In 2SP, the natural way to apply MMT is to integrate it at the tracking stage, replacing the discriminator. 

MMT-integrated DPE would use the code delays estimated from an MMT-aided tracking (`tracking_MMT.m`) to act as the reference code delay for the peak of the ACF. The original DPE_module.m from [DPE_module v1.0](https://github.com/Sergio-Vicenzo/GPSL1-DPEmodule) has been extended to `DPE_module v2.0` perform this feature. 

Since DPE traditionally does not require estimation of code delays, the MMT-integrated DPE is introduced as a variant of DPE instead, one that is specifically designed for urban environments. 

Parts of SoftGNSS v3.0 have been adapted to allow both `tracking_MMT.m` and `DPE_module v2.0` to be integrated into it. This software can run both conventional DPE as well as MMT-integrated DPE. 

Both `tracking_MMT.m` and `DPE_module v2.0` can be integrated into any 2SP MATLAB Software Defined Receivers (SDRs) the same way as [DPE_module v1.0](https://github.com/Sergio-Vicenzo/GPSL1-DPEmodule), as illustrated below.

![DPE plug-in module figure](https://github.com/Sergio-Vicenzo/GPSL1-MMT-DPEmodule/blob/main/DPE%20plug-in%20module%20figure.jpg)

## Running the software

The software presented in this repository is a SoftGNSS v3.0 that has been integrated with `tracking_MMT.m` and `DPE_module v2.0`, and using it follows the same steps as running a regular SoftGNSS v3.0. Modifications were made to SoftGNSS v3.0 to also run on 16-bit samples (`int16` or `short` data types). 

When running with 16-bit data samples, use `int16` in `settings.dataType` (instead of `short`). On the other hand, use `schar` in `settings.dataType` (instead of `int8`) when running with 8-bit data samples.

Further information on using the software can be found in the original SoftGNSS v3.0 readme `readme - SoftGNSS.txt` and ppt `GPS_L1_CA_SDR.pdf`.

## MMT-DPE_module configuration
All `DPE_module v2.0` and `tracking_MMT.m` parameters can be edited from `initSettings.m`, which are listed below.

1. `settings.candPVT_spacing` = Grid spacing for the latitude-longitude-height estimation (Default = 1 meter)

2. `settings.DPE_latlong_span` = Span of latitude and longitude search space (Default = ±30 meters)

3. `settings.DPE_height_span` = Span of height search space (Default = ±50 meters)

4. `settings.DPE_clkBias_span` = Span of clock bias search space (Default = ±20 meters)

5. `settings.DPE_nonCohInt` = DPE non-coherent integration time (Default = 1 ms) to improve the performance of DPE as the satellite correlations would be more filtered

6. `settings.DPE_plotCorrelogram` = Output the correlograms, plotted at the estimated DPE height

7. `settings.gt_llh` = Ground truth coordinates in geodetic coordinates to output the positioning errors from both DPE and 2SP

8. `settings.chipspacing_dpe_precalc` = Chip spacings between the pre-calculated correlations (Default = chips/sample)

9. `settings.MMT` = 1 (Activate MMT-integrated DPE and 2SP) or 0 (Perform conventional DPE and 2SP)

10. `settings.MMT_const ` = Relative amplitude constraint for MMT (Default = 0.8)

## Dependencies

This software requires the `Parallel Computing Toolbox` from MATLAB to accelerate MMT computation at `tracking_MMT.m`.

`tracking_MMT.m` and `DPE_module v2.0` were developed with MATLAB2022a and it is recommended to run the program with the same version or later. Running the program with other MATLAB versions has yet to be tested. No additional MATLAB toolbox is required i.e., it can run with just the basic MATLAB package and no additional toolboxes.

Dependencies of the underlying SoftGNSS v3.0 can be found in its readme file `readme - SoftGNSS.txt`

## Sample IF datasets
Sample GPS L1 C/A intermediate frequency (IF) datasets can be accessed through the following link

[Sample IF datasets](https://drive.google.com/drive/folders/12i75AUCq3DoXqF6xqQ88tibIY3eSlucN?usp=sharing)

All the datasets were collected in a static urban scenario at The Hong Kong Polytechnic University main campus. The collection site taken with Google Earth is shown below (marked H).

![LABSAT data collection site](https://github.com/Sergio-Vicenzo/GPSL1-DPEmodule/blob/main/Collection%20Site.jpg)

The datasets have the following configuration.

`FRONT-END EQUIPMENT`		: LABSAT 3 Wideband

`SAMPLING FREQUENCY`		: 58 MHz				(settings.samplingFreq = 58e6)

`INTERMEDIATE FREQUENCY`	: 1580e6-1575.42e6 Hz 			(settings.IF = 1580e6-1575.42e6)

`DATA TYPE`			: 'schar' or 'int8' 			(settings.dataType = 'schar')

`FILE TYPE`			: IQ 					(settings.fileType = 2)

`FRONT-END BANDWIDTH`		: 56 MHz

`GROUND TRUTH` (in llh)		: [22.30397720 114.17889380 10.770] 	(settings.gt_llh = [22.30397720 114.17889380 10.770])

`FILE LENGTH`			: 240,000 ms 				(settings.msToProcess = 240000)


## Author

`tracking_MMT.m` and `DPE_module v2.0` was written by Sergio Vicenzo.

Sergio Vicenzo is currently a PhD candidate at the Department of Aeronautical and Aviation Engineering, The Hong Kong Polytechnic University. He received first-class honours in Bachelor of Engineering in Aviation Engineering from the same university in 2022. He is supervised by Assistant Professor Bing Xu, and co-supervised by Associate Professor Li-Ta Hsu. His research interests include GNSS urban navigation and positioning with direct position estimation.

Email: <seergio.vicenzo@connect.polyu.hk>

## Disclaimer
I, Sergio Vicenzo, hereby do not claim ownership nor any responsibility for the SoftGNSS v3.0 software by Prof. Kai Borre and Prof. Dennis M. Akos and the `llh2xyz.m` script by Prof. Todd Walter. `tracking_MMT.m` is extended from `tracking.m` from SoftGNSS v3.0.

All changes made to the SoftGNSS v3.0 software (including `tracking_MMT.m`) have been labelled appropriately and its GNU license has been kept intact in each code file, in accordance to GNU General Public License. No changes were made to the `llh2xyz.m` script and its copyright statement is kept intact.

No copyright infringement is intended from any of the MATLAB scripts provided in this GitHub repository.

## References

Tang S, Li H, Closas P (2024) Assessment of Direct Position Estimation Performance in Multipath Channels. In: Proc. ION GNSS+ 2024, Institute of Navigation, Baltimore, Maryland, USA, September 16 - 20, 3705 - 3714.

Vicenzo S, Xu B, Dey A, Hsu L-T (2023) Experimental Investigation of GNSS Direct Position Estimation in Densely Urban Area. In: In: Proc. ION GNSS+ 2023, Institute of Navigation, Denver, Colorado, USA, September 19 – 23, 2906-2919.

Weill LR (2002) Multipath Mitigation using Modernized GPS Signals: How Good Can it Get?. In: Proc. ION GPS 2002, Institute of Navigation, Portland, Oregon, USA, September 24 - 27, 493 - 505.

## Note
This repository is still under development. I apologise if some information is incomplete or if you encounter any problems...

If you have any questions, drop an email at <seergio.vicenzo@connect.polyu.hk>!

Last Updated: 05 Nov 2024

	   
