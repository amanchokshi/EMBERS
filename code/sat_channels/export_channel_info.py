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

    ref_tiles = ['rf0XX', 'rf0YY', 'rf1XX', 'rf1YY']

    for ref_tile in ref_tiles:

        Path(f'{out_dir}/{ref_tile}/histograms').mkdir(parents=True, exist_ok=True)

        norad_ids = []
        channels = []

        ref_dir = Path(f'{out_dir}/{ref_tile}/channel_data')
        
        for f in ref_dir.glob('*.json'):
            with open(f) as chan:
                chan_data = json.load(chan)

                sat_id = chan_data["sat_id"]
                pop_chans = chan_data["pop_chans"]
                sat_count = chan_data["sat_count"]
                chans = chan_data["chans"]
                sats = chan_data["sats"]

                norad_ids.append(sat_id)
                channels.append(pop_chans)

                if chans != []:
            
                    values, counts = np.unique(chans, return_counts=True)

                    #plt_hist(out_dir, values, counts, norad_id)
                    plt_hist(values, counts,
                            'Channel Number',
                            'Number of Passes', 
                            f'Probable Transmission Channel of Sat {sat_id}: {pop_chans}',
                            f'{out_dir}/{ref_tile}/histograms/{sat_id}_channels_histo_{sat_count}_passes_{pop_chans}.png',
                            'GnBu_d')
                    
                if sats != []:
                    s_values, s_counts = np.unique(sats, return_counts=True)
            
                    # only plot sat if it has more than 3 counts
            
                    s_values = s_values[np.where(s_counts>2)]
                    s_counts = s_counts[np.where(s_counts>2)]
                    
                    #plt_hist(out_dir, values, counts, norad_id)
                    plt_hist(s_values, s_counts,
                            'Norad Catalogue ID',
                            'Number of Passes', 
                            f'Possible Satellites in {sat_id} Window',
                            f'{out_dir}/{ref_tile}/histograms/{sat_id}_sats_histo.png',
                            'rocket')


        sat_channels = dict(zip(norad_ids, channels))
                    
        Path(f'{out_dir}/channel_maps').mkdir(parents=True, exist_ok=True)
        
        with open(f'{out_dir}/channel_maps/{ref_tile}_channel_map.json','w') as f: 
            json.dump(sat_channels, f) 

