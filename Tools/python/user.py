import os

if os.environ["USER"] in ["robert.schoefbeck"]:
    #postprocessing_output_directory = "/mnt/hephy/cms/robert.schoefbeck/tWZ/nanoTuples"
    postprocessing_output_directory = "/scratch-cbe/users/robert.schoefbeck/tWZ/nanoTuples"
    postprocessing_tmp_directory    = "/scratch/hephy/cms/robert.schoefbeck/tWZ/tmp/"
    plot_directory                  = "/mnt/hephy/cms/robert.schoefbeck/www/tWZ/plots"
    cache_dir                       = "/mnt/hephy/cms/robert.schoefbeck/tWZ/caches"
    # Analysis result files
    analysis_results                = "/mnt/hephy/cms/robert.schoefbeck/tWZ/results/v1"
    mva_directory                   = "/mnt/hephy/cms/robert.schoefbeck/tWZ/MVA"
    cern_proxy_certificate          = "/users/robert.schoefbeck/.private/.proxy"

if os.environ["USER"] in ["suman.chatterjee"]:
    histogram_output_directory = "/scratch-cbe/users/suman.chatterjee/ParticleNet/Histograms"
    postprocessing_output_directory = "/scratch-cbe/users/suman.chatterjee/VH/Dilep/"
    postprocessing_tmp_directory    = "/scratch/hephy/cms/suman.chatterjee/VH/tmp/"
    #plot_directory                  = "/mnt/hephy/cms/suman.chatterjee/www/VH/plots"
    #cache_dir                       = "/mnt/hephy/cms/suman.chatterjee/Higgs_EFT/caches"
    # Analysis result files
    #analysis_results                = "/mnt/hephy/cms/suman.chatterjee/Higgs_EFT/results"
    #mva_directory                   = "/mnt/hephy/cms/suman.chatterjee/Higgs_EFT/MVA"
    #cern_proxy_certificate          = "/users/suman.chatterjee/.private/.proxy"

if os.environ["USER"] in ["andreas.gruber"]:
   histogram_output_directory = "/scratch-cbe/users/andreas.gruber/ParticleNet/Histograms"
