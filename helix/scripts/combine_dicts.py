import glob, pickle, json, sys, os


if os.path.exists('final.pkl'):
    print('final.pkl exists; opening...')
    with open('final.pkl', 'rb') as f:
        out_dict = pickle.load(f)
else:
    out_dict = {}
if os.path.exists('parsed.json'):
    print('Some dictionaries were previously parsed. Opening list of parsed dictionaries.')
    with open('parsed.json', 'r') as f:
        parsed = json.load(f)

else:
    parsed = []
for path in glob.glob('*.pkl'):
    if path not in parsed and path!='final.pkl':
        try:
            print('Adding {}'.format(path))
            with open(path, 'rb') as f:
                table = pickle.load(f)
            for key in table:
                if key not in out_dict:
                    out_dict[key] = {}
                out_dict[key].update(table[key])
            print('{} is taking up {}G of memory'.format(
                path, sys.getsizeof(table) * (10**-9)
                ))
            del table
            print('Final dict is taking up {}G of memory'.format(
                sys.getsizeof(out_dict) * (10**-9)
                ))
            parsed.append(path)
        except KeyboardInterrupt:
            print("Keyboard interrupt; saving progress...")
            with open('final.pkl', 'wb') as f:
                pickle.dump(out_dict, f)
                del out_dict
            with open('parsed.json', 'w') as f:
                json.dump(parsed, f)

with open('final.pkl', 'wb') as f:
    pickle.dump(out_dict, f)
