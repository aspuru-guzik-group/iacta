import io_utils
import analysis
from analysis import *
import sys
import pandas as pd
import shutil, os
from constants import hartree_ev, ev_kcalmol
import json


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Take a reaction from minigabe, and extract the quantities of interest for further analysis."
        )
    parser.add_argument("reaction_id",
                        help="Reaction # as obtained from read_reactions.py or from "+
                        "looking at reactions.csv.", type=int)
    parser.add_argument("folder", help="Folder of the reaction search run.")
    parser.add_argument("reaction_folder", help="Folder where reaction results are stored.")
    parser.add_argument("-w", help="Overwrite reaction_folder.", action='store_true')

    args = parser.parse_args()
    folder =args.folder
    out_dir = args.reaction_folder


    print("Loading reaction data from %s ..." % folder)
    species = pd.read_csv(folder + "/species.csv").set_index('smiles')
    reactions = pd.read_csv(folder + "/reactions.csv").set_index('reaction_id')
    trajectories = pd.read_csv(folder + "/trajectories.csv")

    print("Setting up output directory %s " % out_dir)
    try:
        os.makedirs(out_dir)
    except FileExistsError:
        print("Output directory exists!")
        if args.w:
            # Delete the directory, make it and restart
            print("   👍 but that's fine! -w flag is on.")
            print("   📁 %s is overwritten."% out_dir)
            shutil.rmtree(out_dir)
            os.makedirs(out_dir)
        else:
            print("   👎 -w flag is off -> exiting! 🚪")
            raise SystemExit(-1)

    reaction = reactions.loc[args.reaction_id]
    reactant = species.loc[reaction['from']]
    product = species.loc[reaction['to']]
    trajs = trajectories[trajectories.reaction_id == args.reaction_id]
    trajs = trajs.drop(columns='reaction_id').set_index('traj_id')

    logf = open(out_dir + "/summary", "w")
    def print(s):
        sys.stdout.write(s + "\n")
        logf.write(s + "\n")


    TS = (trajs.E_TS.min() - reactant.E) * hartree_ev * ev_kcalmol
    dE = (product.E - reactant.E) * hartree_ev * ev_kcalmol
    print("------------------------------------------------------------------------------")
    print("                                Reaction #%i"% args.reaction_id)
    print("from:\n    %s" % reactant.name)
    print("  to:\n    %s" % product.name)
    print("wasΔsampled in %i trajectories." % len(trajs))
    print("\nTrajectory-derived energies:")
    print("    E(R->P)  = %8.4f kcal/mol" % dE)
    print("    E(R->TS) = %8.4f kcal/mol" % TS)

    print("")
    print("    Reactant conformer: %s" % reactant.file )
    print("       -> reactant.xyz")
    shutil.copy(reactant.file, out_dir + "/reactant.xyz")
    print("     Product conformer: %s" % product.file )
    print("       -> product.xyz")
    shutil.copy(product.file, out_dir + "/product.xyz")
    print("\n\n                          Trajectories")
    print("-----+------------------------------------------------------+-----------------")
    print("     |   Energies of important conformers in kcal / mol     | Number of")
    print("  #  |      R        R*        TS        P*         P       | intermediates")

    for k, row in trajs.iterrows():
        fn = row.folder + "/reaction_data.json"
        with open(fn,"r") as fin:
            react = json.load(fin)

        path_data = [(ind,p) for ind,p in enumerate(react['stretch_points'])
                     if ((p >= row.i_R) and (p <= row.i_P))]
        path_indices, cd_points = zip(*path_data)
        start = min(cd_points)
        end = max(cd_points)+ 1

        full_path, full_E = io_utils.traj2str("%s/opt.xyz" % row.folder)
        full_path = full_path[start:end]
        full_E = full_E[start:end]

        if row.i_P < row.i_R:
            # The reaction is actually the reverse of the CD trajectory.
            path_indices = path_indices[::-1]
            cd_points = cd_points[::-1]
            full_path = full_path[::-1]
            full_E = full_E[::-1]

        path_out = "%s/%4.4i" % (out_dir, row.name)
        with open("%s_cd_path.xyz" % path_out, "w") as f:
            for s in full_path:
                f.write(s)

        ninter = 0
        nts =0
        with open("%s_states.xyz" % path_out, "w") as f:
            with open(out_dir + "/reactant.xyz", "r") as struct:
                shutil.copyfileobj(struct, f)

            for i, pi in zip(path_indices, cd_points):
                if react['is_stable'][i]:
                    with open(row.folder + "/stable_%4.4i.xyz"%pi, "r") as struct:
                        shutil.copyfileobj(struct, f)
                    ninter += 1
                else:
                    with open(row.folder + "/ts_%4.4i.xyz"%pi, "r") as struct:
                        shutil.copyfileobj(struct, f)
                    nts += 1

            with open(out_dir + "/product.xyz", "r") as struct:
                shutil.copyfileobj(struct, f)

        shutil.copy(row.folder + "/ts_%4.4i.xyz" % row.i_TS, "%s_TS.xyz" % path_out)


        conv = hartree_ev * ev_kcalmol
        print("%4.4i | %+8.3f  %+8.3f  %+8.3f  %+8.3f  %+8.3f     |    %2i "
              % (row.name,
                 0.0,
                 (full_E[0] - reactant.E) * conv,
                 (row.E_TS - reactant.E) * conv,
                 (full_E[-1] - reactant.E) * conv,
                 (product.E - reactant.E) * conv,
                 ninter-2))


    print("-----+------------------------------------------------------+--------------")
    print("     | R  = Reactants                                       |              ")
    print("     | R* = Activated reactants on the trajectory           |              ")
    print("     | TS = Transition state                                |              ")
    print("     | P* = Activated products on the trajectory            |              ")
    print("     | P  = Products                                        |              ")
    print("     +------------------------------------------------------+              ")
    logf.close()
