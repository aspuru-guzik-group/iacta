import io_utils
import analysis
from analysis import *


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Read reaction products."
        )
    parser.add_argument("folder", help="Folder containing the react*** files.")
    parser.add_argument("-c", "--resolve-chiral",
                        help="Resolve chiral compounds. Defaults to no.",
                        action="store_true")
    parser.add_argument("--all", help="Do not limit reaction network to reactant"
                        +" and connected products. Defaults to false.",
                        action="store_true")
    parser.add_argument("--ts", help="Sort products by TS energy as opposed"
                        +" to enthalpy. (the default).",
                        action="store_true")
    parser.add_argument("--local", help="Use reaction-local barrier instead of TS energy.",
                        action="store_true")
    args = parser.parse_args()
    folder =args.folder
    pathways = read_all_reactions(folder)
    species = get_species_table(pathways, resolve_chiral=args.resolve_chiral)


    reactant, E = io_utils.traj2smiles(folder + "/init_opt.xyz", index=0)
    if args.all:
        final = analyse_reaction_network(pathways,species,list(species.index),
                                         sort_by_barrier=args.ts,
                                         reaction_local=args.local,
                                         resolve_chiral=args.resolve_chiral)
    else:
        print("Reactant: %s" % reactant)
        if reactant in species.index:
            final = analyse_reaction_network(pathways,species,[reactant],
                                             sort_by_barrier=args.ts,
                                             reaction_local=args.local,
                                             resolve_chiral=args.resolve_chiral)
        else:
            print("Error! Reactant not in found species")
            raise SystemExit(-1)

    # Finally, save parsed reaction network
    final.to_csv(folder + "/parsed_reactions.csv")
    species.to_csv(folder + "/parsed_species.csv")
