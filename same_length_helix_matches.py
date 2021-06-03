# same_length_match_indices.py

import os
import pandas as pd
import pickle as pkl

path = 'same_length_matcher_outputs/'
dirlist = os.listdir(path)
print(dirlist)

rows = []
successes = [] # make this a list of lists containing: [pdbid, 3turn or 'no', 4turn or 'no', etc]
for d in dirlist:
	idx = 0
	success_entry = [d.split('_')[0].lower(), 'no_3', 'no_4', 'no_6', 'no_8']
	print(d)
	subdirlist = os.listdir(path + d + '/') # these should be 3, 4, 6, and 8turn
	subdirlist.sort()
	for subdir in subdirlist: # should be of length 4
		idx += 1
		print(subdir)
		# open_file = open(path + d + '/query_results_000.pkl', "rb")
		# curr_df = pkl.load(open_file)
		# open_file.close()
		curr_df = pd.read_pickle(path + d + '/' + subdir + '/query_results_000.pkl')
		query_helices = os.path.join('./rifdock',
				'all_outputs', d, 'cluster_representatives',
				subdir, 'query_helices.pkl')
		curr_df['query_helices'] = query_helices
		print(curr_df)
		row = curr_df.loc[lambda curr_df: curr_df['name'].str.startswith(d.split('_')[0].lower())]
		if not row.empty:
			success_entry[idx] = subdir
			rows.append(row)
	successes.append(success_entry)
	print(success_entry)

	# Add to the pkl file after each directory, in case it freezes again
	open_file = open('same_length_match_indices.p', "wb")
	pkl.dump(rows, open_file)
	open_file.close()

	open_file = open('same_length_successes.p', "wb")
	pkl.dump(successes, open_file)
	open_file.close()

print(rows)

# open_file = open('same_length_match_indices.p', "wb")
# pkl.dump(rows, open_file)
# open_file.close()

# open_file = open('same_length_successes.p', "wb")
# pkl.dump(successes, open_file)
# open_file.close()

# output_df = pd.DataFrame(rows) # review_matches is expecting a df in this format

# output_df.to_pickle('052521match_indices.p')
