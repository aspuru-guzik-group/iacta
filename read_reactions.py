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
    args = parser.parse_args()
    folder =args.folder
    pathways = read_all_reactions(folder, resolve_chiral=args.resolve_chiral)
    species = get_species_table(pathways)


    reactant, E = io_utils.traj2smiles(folder + "/init_opt.xyz", index=0)
    if args.all:
        final = analyse_reaction_network(pathways,species,list(species.index))
    else:
        final = analyse_reaction_network(pathways,species,[reactant])

    # Finally, save parsed reaction network
    final.to_csv(folder + "/parsed_reactions.csv")
    species.to_csv(folder + "/parsed_species.csv")    

