if __name__=="__main__":

    import sys
    import json
    import argparse
    import numpy as np
    from pathlib import Path
    from channels_plt import plt_hist
    
    
    parser = argparse.ArgumentParser(description="""
        Determines which channel each satellite occupies
        """)
    
    parser.add_argument('--chan_dir', metavar='\b', default='../../outputs/sat_channels/channel_data', help='Dir where channel data json files are saved. Default=../../outputs/sat_channels/channel_data')
    parser.add_argument('--out_dir', metavar='\b', default='./../../outputs/sat_channels',help='Output directory. Default=./../../outputs/sat_channels')
    

    
    args = parser.parse_args()
    
    chan_dir =          Path(args.chan_dir)
    out_dir =           Path(args.out_dir)


    Path(f'{out_dir}/histograms').mkdir(parents=True, exist_ok=True)

    chan_dir = Path(f'{out_dir}/channel_data')
    
    for f in chan_dir.glob('*.json'):
        with open(f) as chan:
            chan_data = json.load(chan)

            sat_id = chan_data["sat_id"]
            pop_chans = chan_data["pop_chans"]
            chans = chan_data["chans"]

            if chans != []:
        
                values, counts = np.unique(chans, return_counts=True)

                #plt_hist(out_dir, values, counts, norad_id)
                plt_hist(values, counts,
                        'Channel Number',
                        'Number of Passes', 
                        f'Probable Transmission Channel of Sat {sat_id}: {pop_chans}',
                        f'{out_dir}/histograms/{sat_id}_channels_histo_{sum(counts)}_passes_{pop_chans}.png',
                        'GnBu_d')
                

