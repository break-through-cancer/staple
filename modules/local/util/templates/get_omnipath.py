#!/usr/bin/env python
from omnipath.interactions import import_intercell_network
import pandas as pd
import logging as log
import json
from importlib.metadata import version
log.basicConfig(level=log.INFO)

if __name__ == "__main__":
    interactions_params = json.loads('${params.om_interaction_params}') if '${params.om_interaction_params}' else None
    receiver_params = json.loads('${params.om_receiver_params}') if '${params.om_receiver_params}' else None
    transmitter_params = json.loads('${params.om_transmitter_params}') if '${params.om_transmitter_params}' else None
    interaction_query = '${params.om_interaction_query}' if '${params.om_interaction_query}' else None

    # Fetch the default intercell/ligand-receptor network used by Squidpy
    log.info("Fetching Omnipath intercell network...")
    interactions = import_intercell_network(interactions_params=interactions_params,
                                            receiver_params=receiver_params,
                                            transmitter_params=transmitter_params)
    # omnipath returns protein names, rename
    log.info(f"Setting source and target to gene symbols...")
    interactions.rename(columns={"source": "omnipath_source", "target": "omnipath_target"}, inplace=True)
    interactions.rename(columns={"genesymbol_intercell_source": "source", "genesymbol_intercell_target": "target"}, inplace=True)
    
    # repeat what squidpy does with complex names
    interactions["source"] = interactions["source"].str.replace("^COMPLEX:", "", regex=True)
    interactions["target"] = interactions["target"].str.replace("^COMPLEX:", "", regex=True)

    # sort to keep source and target as first two columns
    interactions = interactions[["source", "target"] + [col for col in interactions.columns if col not in ["source", "target"]]]
    
    #apply additional filters if requested
    if interaction_query:
        try:
            filtered_interactions = interactions.query(interaction_query)
            if filtered_interactions.empty:
                log.warning(
                    f"No interactions left after applying query '{interaction_query}'; "
                    "keeping unfiltered interactions."
                )
            else:
                interactions = filtered_interactions
        except Exception as e:
            log.warning(
                f"Failed to apply custom interaction query '{interaction_query}': {e}; "
                "keeping unfiltered interactions."
            )

    # Save to CSV
    log.info(f"Saving {len(interactions)} interactions to omnipath.csv...")
    interactions.to_csv("omnipath.csv", index=False)
    log.info("Omnipath intercell network saved successfully.")
    
    #versions
    with open("versions.yml", "w") as f:
        f.write("${task.process}:\\n")
        f.write("    omnipath: {}\\n".format(version("omnipath")))