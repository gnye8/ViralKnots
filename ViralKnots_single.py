from ViralKnots_utils import *

if __name__=='__main__':

    parser=argparse.ArgumentParser()

    parser.add_argument('--seqs', nargs='+', required=True)
    parser.add_argument('--coords', nargs='+', required=True)
    parser.add_argument('--window', '-w', type=int, required=True)
    parser.add_argument('--pk_predictors', nargs='+', required=True)
    parser.add_argument('--bpp_packages', default=[], nargs='+')
    parser.add_argument('--temp_folder', default=None, type=str)
    parser.add_argument('--linear_partition', default=True, type=bool)
    args=parser.parse_args()

    struct_list = []
    for pk_predictor in args.pk_predictors:
        if pk_predictor == 'threshknot':
            for bpp_package in args.bpp_packages:
                for seq,coord in zip(args.seqs, args.coords):
                    struct_list.append(get_structure(seq, coord, pk_predictor, args.window, bpp_package, args.linear_partition))
        else:
            for seq,coord in zip(args.seqs, args.coords):
                struct_list.append(get_structure(seq, coord, pk_predictor, args.window, bpp_package=None, linear_partition=args.linear_partition))

    df = pd.DataFrame(struct_list,columns=["predictor","start","end","sequence", "struct", "pseudoknot"])
    #TO DO: fix the way these are named so that you can use a smaller job size without it erroring because the file name is too long
    csv_name = ''
    for i in range(len(args.coords)):
        csv_name += args.coords[i] + '_'

    df.to_csv(f"{args.temp_folder}/{csv_name}.csv", index=False)
