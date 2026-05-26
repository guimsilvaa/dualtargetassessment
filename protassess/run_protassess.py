import os
import logging
from datetime import datetime
from detect_cavities import (process_target_folder)
from extract_descriptors import (build_descriptor_dataframe)
from analyze_ensembles import (analyze_ensembles)
from cluster_pockets import (cluster_pockets)

def main():
    # Detect target folders
    target_folders = sorted([
        f for f in os.listdir(".")
        if (
            os.path.isdir(f)
            and f.startswith("target")
        )
    ])
    if len(target_folders) < 2:
        raise ValueError(
            "At least two target folders "
            "are required."
        )

    # Output directory
    timestamp = datetime.now().strftime(
        "%y%m%d"
    )
    output_dir = (
        f"{timestamp}_"
        f"{'_'.join(target_folders)}"
        f"_protassess-out"
    )
    os.makedirs(
        output_dir,
        exist_ok=True
    )

    # Logging
    log_file = os.path.join(
        output_dir,
        "protassess.log"
    )
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format=(
            "%(asctime)s - "
            "%(levelname)s - "
            "%(message)s"
        )
    )

    # Header
    print("\n===================================")
    print(" ProtAssess ")
    print(" Ensemble Pocket Analysis Tool ")
    print("===================================\n")
    print("Detected targets:")
    for t in target_folders:
        print(f" - {t}")
    print(
        f"\nOutput directory:\n{output_dir}\n"
    )
    logging.info(
        "ProtAssess started."
    )
    logging.info(
        f"Targets: {target_folders}"
    )

    # Process all targets
    all_results = []
    for target_folder in target_folders:
        results = process_target_folder(
            target_folder,
            output_dir
        )
        all_results.extend(results)

    # Build descriptor dataframe
    descriptor_df = (
        build_descriptor_dataframe(
            all_results
        )
    )
    descriptors_csv = os.path.join(
        output_dir,
        "descriptors.csv"
    )
    descriptor_df.to_csv(
        descriptors_csv,
        index=False
    )
    logging.info(
        f"Descriptor CSV saved: "
        f"{descriptors_csv}"
    )
    print(
        "Descriptor extraction complete."
    )

    # Ensemble analysis
    analyze_ensembles(
        descriptors_csv,
        output_dir
    )
    logging.info(
        "Ensemble analysis complete."
    )
    print(
        "Ensemble analysis complete."
    )

    # Clustering
    cluster_pockets(
        descriptors_csv,
        output_dir
    )
    logging.info(
        "Clustering complete."
    )
    print(
        "Pocket clustering complete."
    )

    # Final message
    logging.info(
        "ProtAssess finished successfully."
    )
    print("\n===================================")
    print(" ProtAssess Finished Successfully ")
    print("===================================\n")
    print(
        f"Results saved in:\n{output_dir}\n"
    )
if __name__ == "__main__":
    main()