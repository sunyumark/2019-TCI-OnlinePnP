# [An Online Plug-and-Play Algorithm for Regularized Image Reconstruction](https://ieeexplore.ieee.org/document/8616843/)

The arxiv version of this paper is available [here](https://arxiv.org/abs/1809.04693)

## Abstract
Plug-and-play priors (PnP) is a powerful framework for regularizing imaging inverse problems by using advanced denoisers within an iterative algorithm. Recent experimental evidence suggests that PnP algorithms achieve the state-of-the-art performance in a range of imaging applications. In this paper, we introduce a new online PnP algorithm based on the proximal gradient method (PGM). The proposed algorithm uses only a subset of measurements at every iteration, which makes it scalable to very large datasets. We present a new theoretical convergence analysis, for both batch and online variants of PnP-PGM, for denoisers that do not necessarily correspond to proximal operators. We also present simulations illustrating the applicability of the algorithm to image reconstruction in diffraction tomography. The results in this paper have the potential to expand the applicability of the PnP framework to very large datasets.

## PnP-OPGM vs. PnP-PGM vs. PnP-ADMM for diffraction tomogaphy
![visualExamples](images/compare.gif)

## Citation
    @ARTICLE{8616843,
    author={Y. {Sun} and B. {Wohlberg} and U. S. {Kamilov}},
    journal={IEEE Transactions on Computational Imaging},
    title={An Online Plug-and-Play Algorithm for Regularized Image Reconstruction},
    year={2019},
    volume={5},
    number={3},
    pages={395-408},
    doi={10.1109/TCI.2019.2893568},
    ISSN={2573-0436},
    month={Sep.},}
