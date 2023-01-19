from ViralKnots_utils import *

if __name__=='__main__':

    parser=argparse.ArgumentParser()

    parser.add_argument('--seqs', nargs='+', required=True)
    parser.add_argument('--coords', nargs='+', required=True)
    parser.add_argument('--shape_windows', nargs='+', required=True)
    parser.add_argument('--window', '-w', type=int, required=True)
    parser.add_argument('--shape_name', type=str, required=True)
    parser.add_argument('--temp_folder', default=None, type=str)

    args=parser.parse_args()

    shape_windows_str = args.shape_windows
    shape_windows = []
    for shape_window_str in shape_windows_str:
        shape_windows.append(list(map(float, shape_window_str.split(','))))

    struct_list = []
    for seq,coord,shape in zip(args.seqs,args.coords,shape_windows):
        struct_list.extend(shapeknots_predict(seq, shape, args.shape_name, coord, args.window))

    df = pd.DataFrame(struct_list,columns=["predictor","start","end","sequence", "struct", "pseudoknot"])

    csv_name = ''
    for i in range(len(args.coords)):
        # RCK added so does not overwrite pk_predict outputs
        csv_name += args.coords[i] + '_'

    df.to_csv(f"{args.temp_folder}/{csv_name}_{args.shape_name}_shapeknots.csv", index=False)
