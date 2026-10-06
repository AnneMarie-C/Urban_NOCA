# Master's Research Project - Winter Survival and Urban Green Space use of the Northern Cardinal

Analysis code for both chapters of an MSc thesis on winter survival and urban green space use by the Northern Cardinal (*Cardinalis cardinalis*) in Montreal, Quebec, Canada.

> Cousineau, A.-M. (2026). *Winter survival and urban green space use of the northern cardinal.* M.Sc. thesis, Department of Natural Resource Sciences, McGill University. Supervisors: K. Elliott & B. Frei.

## Overview
### Chapter 1
Northern cardinals were radio-tracked across three sites spanning an urbanization gradient in Montreal to test how habitat selection scale and habitat associations shift with urbanization intensity. This code:

- Estimates the scale of effect of habitat variables on cardinal space use using Resource Selection Functions, comparing model fit across a range of buffer radius
- Quantifies home range size per individual/site and compares across the urbanization gradient
- Models habitat associations (vegetation vertical structure, paved road density, proximity to bird feeders) against the best-supported scale at each site

**Key findings this code supports:** cardinals selected habitat at a smaller scale (25 m radius) at the most urbanized site versus a larger scale (50 m radius) at the least urbanized site; home ranges were roughly 3x larger at the most urbanized site; and the habitat variables driving space use differed by site urbanization level.

### Chapter 2
Some of the individuals radio-tracked in chapter 1 remained close enough to the receiver stations to be able to investigate their daily activity with regards to winter conditions. This code:

- Determines the onset and end time of activity for every day of detection for every qualifying individuals
- Uses winter weather data to model the effect of a various set of conditions on the timing of cardinals' activity during the wintertime

**Key findings this code supports:** wintering northern cardinals are highly flexible in their daily schedules, delaying activity under persistent cold but paradoxically starting earlier following acute cold snaps, while only extreme cold causes them to shorten their foraging period at the end of the day.

## Repository structure

```
MSc_Project/
├── NOCA_Winter_Habitat_Association/
├── NOCA_Winter_Daily_Activity/
├── Publications/
├── Thesis/
└── README.md
```

## Requirements

- R (version ≥ 4.x)
- Key packages: `mgcv` (GAMs), `sf` (spatial data), `adehabitatHR` or similar (home range estimation), `tidyverse`

## Usage

Currently under construction!

## Data availability

Currently under construction!

## Citation

If you use this code, please cite the thesis above.

## License

This repository is released under [CC0 1.0](LICENSE) — public domain dedication.

## Contact

Visit my GitHub profile for all the links to my socials.
