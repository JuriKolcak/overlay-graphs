# Overlay graphs (OGs)

Formal spec of overlay graphs (OGs).

The python scripts in this repository can be used to convert
the compatible enzymatic reaction mechanisms obtained from the [M-CSA database](https://www.ebi.ac.uk/thornton-srv/m-csa/) into OGs and substrate rules.

## Dependencies
The scripts in this package rely heavily on a graph transformation framework for chemistry called [MØD](https://cheminf.imada.sdu.dk/mod/).

Other dependencies include [NetworkX](https://networkx.org/) and [RDKit](https://rdkit.org/).

All dependencies can be easily installed using conda. You can either install the included `environment.yml` conda environment using the following command:
```commandline
conda env create -f environment.yml
```

Or create your own:
```commandline
conda create -n [your environment name] -c jakobandersen -c conda-forge mod networkx rdkit
```

## Converting mechanisms to OGs

Run the following command to generate OGs for the mechanism listed in a `mechanisms.json` file:
```commandline
mod -f ogs.py
```

The script produces both an `overlay_graphs.json` file with all the overlay graphs, as well as an visual output `summary/summary.pdf` with at most 10 overlay graphs per reaction sequence.
Please note that the production of the visual summary can be relatively slow in case hundreds of mechanisms are being processed. Allowing the process to run in parallel by adding `-j [number of threads]` option may help.

Finally, the input file `mechanisms.json` is expected to contain full molecule GML rules, corresponding to the elementary steps of mechanisms as depicted in [M-CSA](https://www.ebi.ac.uk/thornton-srv/m-csa/).
The [M-CSA](https://www.ebi.ac.uk/thornton-srv/m-csa/) data can be obtained in the appropriate format using our [scraping and conversion tool](https://github.com/chrstf/mcsa_rule_converter.git).

## Converting OGs to substrate rules

Run the following command to generate substrate rules from the overlay graphs in a `overlay_graphs.json` file:
```commandline
mod -f substrate.py
```

The script produces both an `substrate_rules.json` file with all the substrate rules, as well as an visual output `summary/summary.pdf` with all the substrate rules.
Please note that the production of the visual summary can be relatively slow in case hundreds of overlay graphs are being processed. Allowing the process to run in parallel by adding `-j [number of threads]` option may help.
