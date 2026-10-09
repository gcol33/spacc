"""Fetch the papers behind the feature issues into refs/.

Run from anywhere: python refs/retrieval.py
Codes: B = beta diversity, P = phylogenetic, F = functional, N = null models,
E = extrapolation / richness estimation, C = coverage, I = interaction
networks, R = spatial rarefaction software.

Papers from 2024 onwards are not on Sci-Hub; those come from open-access
repository copies (direct_urls) or, for Elsevier, the univie SSO route.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path.home() / ".python"))
from retrieve_papers import retrieve_all

OUTPUT_DIR = Path(__file__).parent

retrieve_all(
    output_dir=OUTPUT_DIR,
    papers=[
        ("B1_Chao_2023_iNEXT_beta3D_coverage_standardised_beta.pdf", "10.1002/ecm.1588"),
        ("C2_Chiu_2023_Sample_coverage_sample_based_abundance.pdf", "10.1002/ecy.4099"),
        ("R2_Alberdi_2019_Hill_numbers_guide_redundancy.pdf", "10.1111/1755-0998.13014"),
        ("R1_Thouverai_2026_Rarefy_spatially_explicit_rarefaction.pdf", "10.1016/j.ecoinf.2026.103803"),
    ],
    direct_urls=[
        # Europe PMC render of the OUP open-access article (PMC11129585)
        ("P1_Rosindell_2024_EvoHeritage.pdf",
         "https://europepmc.org/articles/PMC11129585?pdf=render"),
        # University of Jyvaskyla repository copy of the CC BY-NC article
        ("N1_Ridder_2024_Shift_and_rotate_null_models.pdf",
         "https://jyx.jyu.fi/bitstreams/6bc9d041-1611-4c1e-a7bf-35707c4a1ec8/download"),
        # University of the Azores repository copy
        ("F2_Dias_2025_Functional_diversity_area_relationships.pdf",
         "https://repositorio.uac.pt/bitstreams/a576e698-3e50-4b69-bce4-affd6d1def60/download"),
        # EcoEvoRxiv preprint of the Ecological Monographs paper
        ("F3_Palacio_2025_Intraspecific_trait_variability_FD.pdf",
         "https://ecoevorxiv.org/repository/object/6455/download/12521"),
        # Oxford Research Archive accepted manuscript
        ("E1_Shestopaloff_2025_Accumulation_rate_curve_estimator.pdf",
         "https://ora.ox.ac.uk/objects/uuid:dc5fbbce-a3a2-435f-a3b2-7b2b70373512/files/rb5644t284"),
        # SLU Epsilon repository copy
        ("B2_Jones_2025_Beta_partition_at_two_resolutions.pdf",
         "https://pub.epsilon.slu.se/id/document/20455101"),
        # Publisher open-access PDFs (may sit behind a bot wall)
        ("F1_Wojcik_2024_Aggregated_functional_diversity_K.pdf",
         "https://besjournals.onlinelibrary.wiley.com/doi/pdfdirect/10.1111/2041-210X.14470"),
        ("E2_Zou_2023_Total_expected_species_rarefaction_fit.pdf",
         "https://esajournals.onlinelibrary.wiley.com/doi/pdfdirect/10.1002/ecs2.4363"),
        ("C1_Chiu_2023_Richness_estimator_sampling_without_replacement.pdf",
         "https://besjournals.onlinelibrary.wiley.com/doi/pdfdirect/10.1111/2041-210X.14146"),
        ("I1_Chiu_2023_Network_diversity_incomplete_sampling.pdf",
         "https://royalsocietypublishing.org/doi/pdf/10.1098/rstb.2022.0183"),
    ],
)
