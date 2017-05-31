#' NCI QIN Head and Neck Cancer PET Segmentation Challenge
#' 
#' A dataset containing tumor volumes derived from different segmentation
#' methods applied to PET imaging scans of cacerous head and neck lesions by
#' institutions participating in a challenge conducted by the NCI Quantitative
#' Imaging Network (QIN).
#' 
#' \tabular{rlll}{
#' [,1] \tab method   \tab factor  \tab unique identifier for the segmentation method. \cr
#' [,2] \tab lesion   \tab factor  \tab identifier for the segmented lesion. \cr
#' [,3] \tab operator \tab factor  \tab identifier for the operator performing the segmentation. \cr
#' [,4] \tab volume   \tab numeric \tab segmented volume (ml). \cr
#' }
#' 
#' @source Beichel RR, Smith BJ, Bauer C, Ulrich EJ, Ahmadvand P, Budzevich MM,
#' Gillies RJ, Goldgof D, Grkovski M, Hamarneh G, Huang Q, Kinahan PE,
#' Laymon CM, Mountz JM, Muzi JP, Muzi M, Nehmeh S, Oborski MJ, Tan Y, Zhao B,
#' Sunderland JJ, Buatti JM. \dQuote{Multi-site quality and variability analysis
#' of 3D FDG PET segmentations based on phantom and clinical image data.}
#' \emph{Medical Physics} 2017; 44(2): 479--496.
"hnc"
