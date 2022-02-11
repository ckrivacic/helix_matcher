'''
Utility functions for plotting.
'''
from matplotlib.lines import Line2D
import os, sys
import time
import subprocess
import pymol


class ClickablePlot(object):
    '''
    Object to allow one to click on points to open structures.
    Usage:
    Press Shift+A to begin selecting multiple models.
    Press Shift+C to stop selecting and open a pymol session.
    Click a point when not in selection mode to open just that model.
    '''
    def __init__(self, ax, df, args, workspace, pvp=False):
        self.ax = ax 
        self.df = df
        self.args = args
        self.workspace = workspace
        self.ax.figure.canvas.mpl_connect("pick_event", self._onpick)
        self.ax.figure.canvas.mpl_connect("key_press_event",
                self._on_key_press)
        self.select = False
        self.selected = []
        self.listening = False
        self.current_ind = None
        self.selected_idx = None
        self.pvp = pvp

    def _on_key_press(self, event):
        if event.key == 'A':
            self.select = True
            print('Selection mode active')
        if event.key == 'C':
            self.select = False
            print('Selection mode disabled')
            if len(self.selected) > 0:
                self.open_selected()
            self.selected = []

        if self.listening:
            if event.key == 'enter':
                self.listening = False
                if self.pvp:
                    design_file_x =\
                            self.df.iloc[self.current_ind[int(self.selected_idx)]]['design_file_{}'.format(self.args['--xaxis'])]
                    design_file_x =\
                            self.df.iloc[self.current_ind[int(self.selected_idx)]]['design_file_{}'.format(self.args['--yaxis'])]
                    fpath_x = self.fetch_file(design_file_x)
                    fpath_y = self.fetch_file(design_file_y)
                    self.selected.extend([fpath_x, fpath_y])
                else:
                    design_file = self.df.iloc[self.current_ind[int(self.selected_idx)]]['design_file']
                    fpath = self.fetch_file(design_file)
                    print('\nSelected file at {}'.format(fpath))
                    if not self.select:
                        self.selected = [fpath]
                        self.open_selected()
                    else:
                        self.selected.append(fpath)
            elif event.key == 'backspace':
                self.selected_idx = self.selected_idx[:-1]
                print('\b', end='', flush=True)
            else:
                print(event.key, end='', flush=True)
                self.selected_idx += event.key

    def open_selected(self):
        print('Launching:')
        print(self.selected)
        pymol.cmd.reinitialize()
        interface_path = '/Users/codykrivacic/software/anaconda3/envs/proteindesign/lib/python3.7/site-packages/pmg_tk/startup/interfaceFinder.py'
        sys.path.insert(0, os.path.dirname(interface_path))
        import interfaceFinder
        pymol.cmd.run(interface_path)
        for design in self.selected:
            design_name = os.path.basename(design).split('.')[0]
            pymol.cmd.load(design, design_name)
            interfaceFinder.interfaceResidues(design_name,
                    selName=design_name + '_interface')
            pymol.cmd.show('sticks', design_name + '_interface')
        pymol.cmd.hide('sticks', 'hydro')
        pymol.cmd.color('orange', 'chain B and name c*')
        pymol.cmd.deselect()
        self.selected = []

        pymol.cmd.save('temp.pse')
        subprocess.call(['pymol', 'temp.pse'])
        return

    def _onpick(self, event):
        ind = event.ind

        if len(ind) > 1:
            print('Open which of the following?')
            printstr = ''
            i = 0
            for idx in ind:
                printstr += '{}: '.format(i)
                if self.pvp:
                    printstr += os.path.basename(self.df.iloc[idx]['design_file_{}'.format(self.args['--xaxis'])])
                else:
                    printstr += os.path.basename(self.df.iloc[idx]['design_file'])
                    printstr += ', score: {}\n'.format(
                            self.df.iloc[idx][self.args['--hue']]
                            )
                i += 1
            print(printstr)
            # inp = int(input(printstr))
            self.selected_idx = ''
            self.listening = True
            self.current_ind = ind
            return
            # while self.listening:
                # Not very elegant but it works...
                # time.sleep(0.1)
            # inp = int(self.selected_idx)
            # patchman_file = self.df.iloc[ind[inp]]['design_file']
        else:
            if self.pvp:
                design_file_x =\
                        self.df.iloc[ind[0]]['design_file_{}'.format(self.args['--xaxis'])]
                print(self.df.iloc[ind[0]]['helix_seq_{}'.format(self.args['--xaxis'])])
                print(self.df.iloc[ind[0]]['helix_seq_{}'.format(self.args['--yaxis'])])
                design_file_y =\
                        self.df.iloc[ind[0]]['design_file_{}'.format(self.args['--yaxis'])]
                fpath_x = self.fetch_file(design_file_x)
                fpath_y = self.fetch_file(design_file_y)
                if not self.select:
                    self.selected = [fpath_x, fpath_y]
                    self.open_selected()
                else:
                    self.selected.extend([fpath_x, fpath_y])
            else:
                patchman_file = self.df.iloc[ind[0]]['design_file']
                fpath = self.fetch_file(patchman_file)
                print('Selected file at {}'.format(fpath))
                if not self.select:
                    self.selected = [fpath]
                    self.open_selected()
                else:
                    self.selected.append(fpath)

    def fetch_file(self, patchman_file):
        '''
        Fetch a file (if not found)
        '''
        pathlist = patchman_file.split('/')
        if os.path.basename(self.workspace.root_dir) not in pathlist:
            fpath = os.path.join(
                    self.workspace.root_dir,
                    patchman_file
                    )
        else:
            start = pathlist.index(
                    os.path.basename(self.workspace.root_dir)
                    ) + 1
            fpath = os.path.join(
                    self.workspace.root_dir,
                    *pathlist[start:]
                    )
        if not os.path.exists(fpath):
            if os.path.basename(self.workspace.root_dir) not in pathlist:
                remote_path = os.path.join(
                        os.path.basename(self.workspace.root_dir),
                        patchman_file
                        )
                remote_path = self.workspace.rsync_url + remote_path
            else:
                rsync_url = self.workspace.rsync_url.split(':')[0] + ':'
                remote_path = rsync_url + patchman_file
            cmd = ['scp', remote_path, fpath]
            print('Running command:')
            print(cmd)
            subprocess.call(cmd)
        return fpath
