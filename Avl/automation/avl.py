""" Avl python module 
cd Avl
python3 automation/avl.py <PLANE_NAME> <TEMPLATE_PATH>
"""
import os
import sys
import subprocess
import warnings

# Simple class for the aircraft object in the AVL stuff
class Aircraft():
    """ Aircraft Class """
    def __init__(self, wing, tail, mass):
        """ Build an aircraft with rectangular wing and tail """
        self.wing=wing
        self.tail=tail
        self.mass=mass

# want a function that we give some general aircraft dimensions and a ref point
# and we get a file into our Avl/runs, for now just run from /workspace
class Avl():
    """  """
    def __init__(self, name:str,  aircraft:Aircraft, avl_path=None, avl_run_path='runs'):
        """ create avl wrapper object 
        Arguments:
            name: name of aircraft associated w/ instance (str)
            aircraft: aircraft object (Aircraft)
            avl_path: path to avl root, /Avl (path_like)
            avl_run_path: path to case runs, /runs (path_like)
        """
        self._name=name
        self._aircraft=aircraft
        self._file=os.path.join(avl_run_path, name+'.avl')
        if not avl_path:
            self._avl_path=''
        else:
            self._avl_path=avl_path
        self._avl_runs_path=avl_run_path

    def _write_surface(self, name, nspan, sspace, dainc, ydup=None):
        """ write surface """
        with open(f'{self._file}', 'a') as avl_file:
            avl_file.write('SURFACE                      | (keyword) \n' \
                           f'{name} \n' \
                           '#Nchord    Cspace   [ Nspan Sspace ] \n' \
                           f'{nspan:0.1f}        {sspace:0.1f} \n\n')
            if ydup is not None:
                avl_file.write('YDUPLICATE \n' \
                           f'{ydup:.1f} \n\n')
            # continue writing to file
            avl_file.write('SCALE \n' \
                           '1.0  1.0  1.0 \n\n' \
                           'TRANSLATE \n' \
                           '0.0  0.0  0.0 \n\n' \
                           'ANGLE \n' \
                           f'  {dainc:.3f}                        | dAinc \n')

    def _write_section(self,
                       x,
                       y,
                       z,
                       chord,
                       ainc,
                       nspan,
                       sspace,
                       airfoil,
                       xhinge=None,
                       ctrl_surf_type=None):
        """ write section into avl file """
        with open(f'{self._file}','a') as avl_file:
            avl_file.write('SECTION                                                     |  (keyword) \n' \
                           f'   {x:.3f}    {y:.3f}    {z:.3f}    {chord:.3f}   {ainc:.3f}    {nspan:.0f}    {sspace:.0f}   ' \
                           '| Xle Yle Zle   Chord Ainc   [ Nspan Sspace ] \n\n' \
                           'AFIL 0.0 1.0 \n' \
                           f'{airfoil}\n')
        # write in control surface section
        match ctrl_surf_type:
                # None: no control surface
                case None:
                  return
                # aileron control surface
                case 'aileron':
                    with open(f'{self._file}','a') as avl_file:
                        avl_file.write("CONTROL \n")
                        avl_file.write(f"aileron  1.0  {xhinge:.3f}  0.000  0.000  0.000  -1.0 \n")
                        avl_file.close()
                # elevator control surface
                case 'elevator':
                    with open(f'{self._file}','a') as avl_file:
                        avl_file.write("CONTROL \n")
                        avl_file.write(f"elevator  1.0  {xhinge:.3f}  0.000  0.000  0.000  -1.0 \n")
                        avl_file.close()
                # rudder control surface
                case 'rudder':
                    with open(f'{self._file}','a') as avl_file:
                        avl_file.write("CONTROL \n")
                        avl_file.write(f"rudder  1.0  {xhinge:.3f}  0.000  0.000  0.000  -1.0 \n")
                        avl_file.close()

    def create_avl_file(self):
        """ 
        Create an avl file for specified aircraft
            eventually need a way to intelligently get mesh values

        "'" header """
        # check if the .avl file exists first, then overwrite it
        if os.path.exists(self._file):
            print('file exists, overwriting file')
        
        with open(f'{self._file}', 'w') as avl_file:
            avl_file.write('# AVL File \n')
    
        sref=self._aircraft.wing['span']*self._aircraft.wing['chord']
        cref=self._aircraft.wing['chord']
        bref=self._aircraft.wing['span']
        xref=self._aircraft.mass['x_ref']
        yref=self._aircraft.mass['y_ref']
        zref=self._aircraft.mass['z_ref']
        with open(f'{self._file}','a') as avl_file:
            avl_file.write(f'{self._name} \n'
                           '0.0                                 | Mach \n' \
                           '0     0     0.0                     | iYsym  iZsym  Zsym \n' \
                           f'  {sref:.3f}     {cref:.3f}     {bref:.3f}   | Sref   Cref   Bref \n' \
                           f'  {xref:.3f}     {yref:.3f}    {zref:.3f}   | Xref   Yref   Zref \n' \
                           ' 0.00                               | CDp  (optional) \n \n')
        """ wing section """
        dainc=self._aircraft.wing['incidence']
        self._write_surface('wing', 15.0, 1.0, dainc, 0.0)
        # root of wing
        x=0.0
        y=0.0
        z=0.0
        chord=self._aircraft.wing['chord']
        ainc=0.0
        nspan=15.0
        sspace=3.0
        xhinge=self._aircraft.wing['chord']-self._aircraft.wing['aileron_chord']
        wing_foil=self._aircraft.wing['airfoil']  # be careful with overwriting strings/chars
        self._write_section(x, y, z, chord, ainc, nspan, sspace, wing_foil)
        # create wing at aileron
        y=self._aircraft.wing['span']/2 - self._aircraft.wing['aileron_span']
        nspan=15.0
        sspace=3.0
        self._write_section(x, y, z, chord, ainc, nspan, sspace, wing_foil, xhinge,'aileron')
        # create wing tip
        y=self._aircraft.wing['span']/2
        nspan=15.0
        sspace=3.0
        self._write_section(x, y, z, chord, ainc, nspan, sspace, wing_foil, xhinge, 'aileron')
        
        """ horizontal tail """
        dainc=self._aircraft.tail['incidence']
        self._write_surface('elevator', 10.0, 1.0, dainc, 0.0)
        # root of wing
        x=self._aircraft.tail['l_tail']
        y=0.0
        z=self._aircraft.tail['h_vertical_offset']
        chord=self._aircraft.tail['h_chord']
        ainc=0.0
        nspan=8.0
        sspace=3.0
        xhinge=self._aircraft.tail['h_chord']-self._aircraft.tail['elev_chord']
        h_foil=self._aircraft.tail['h_foil']
        self._write_section(x, y, z, chord, ainc, nspan, sspace, h_foil, xhinge, 'elevator')
        # wing at tip
        y=self._aircraft.tail['h_span']/2
        self._write_section(x, y, z, chord, ainc, nspan, sspace, h_foil, xhinge, 'elevator')

        """ vertical tail/fin """
        self._write_surface('fin', 7.0, 1.0, 0.0)
        # top of fin
        x=self._aircraft.tail['l_tail']
        y=0.0
        z=self._aircraft.tail['v_span']
        chord=self._aircraft.tail['v_chord']
        ainc=0.0
        nspan=8.0
        sspace=3.0
        xhinge=self._aircraft.tail['v_chord']-self._aircraft.tail['rudd_chord']
        v_foil=self._aircraft.tail['v_foil']
        self._write_section(x, y, z, chord, ainc, nspan, sspace, v_foil, xhinge, 'rudder')
        # bottom of fin
        z=0.0
        self._write_section(x, y, z, chord, ainc, nspan, sspace, v_foil, xhinge, 'rudder')

    def run_avl(self, template_path):
        """ load and run the avl files steps 
        this should return or save something useful 
        Arguments: 
            template_path: treating Avl as the working directory (path_like)
        """
        # TODO: #4 change the automatic name to sort of math the template name
        # TODO: #5 use a cleaner place holder in place PLANES in the template 
        steps_path=f'automation/{self._name}_steps.txt'
        self._generate_steps(
            template_path, steps_path, self._name)

        exe=os.path.dirname(self._avl_path)
        with open(f'{steps_path}', 'r') as steps:
            subprocess.run(
                [str(exe)],
                cwd=exe,
                stdin=steps,
                text=True,
                check=True
            )

    @staticmethod
    def _generate_steps(template_path, output_path, plane_name):
        """ generate steps from a template 
        Arguments:
            template_path: path of template (path_like)
            output_path: path of output (path_like)
            pane_name (str/char)
        """
        # check if the template exists
        if os.path.exists(template_path):
            warnings.warn(
                f'template not found: {template_path}',
                category=UserWarning,
                stacklevel=2
            )

        with open(template_path, 'r') as template:
            # read steps from template
            steps=template.read()
        # replace key with plane_name
        steps=steps.replace('PLANE', plane_name)
        # make a new file where we want the steps
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        # check if the new output exists
        if os.path.exists(output_path):
            print(f'overwriting file at {output_path} \n')

        # write the new steps into the new file
        with open(output_path, 'w') as steps_out:
            steps_out.write(steps)

    # frame: xyz=BLU
# use a dictionary since we can use a yaml for this and we can just generate a ton of yamls
# use avl/runs/yaml?
wing={
    'span': 1.524,
    'chord': 0.416,
    'airfoil': 'clark_y.dat',
    'aileron_span': 0.381,   # span of EACH aileron
    'aileron_chord': 0.104,
    'incidence': 2.0
}

tail={
    'l_tail': 1.5,
    'h_span': 0.5467,
    'h_chord': 0.1647,
    'h_foil': 'naca_0012.dat',
    'h_vertical_offset': 0.1127,
    'elev_chord': 0.08,
    'v_span': 0.2254,
    'v_chord': 0.1543,
    'v_foil': 'naca_0012.dat',
    'v_vertical_offset': None,  # if the tail is offset (skipping for now)
    'rudd_chord': 0.08,
    'incidence': -2.0 # i think idk
}

mass={  # cog of our plane, well have multiple cases, so we can make case yamls
    'x_ref': 0.1543,
    'y_ref': 0.0,
    'z_ref': 0.0,
}

# aircraft=Aircraft(wing, tail, mass)
# avl=Avl(name='test_plane', aircraft=aircraft)
# avl.create_avl_file()

def main(*args):
    """ thing to run at call """
    aircraft=Aircraft(wing, tail, mass)
    avl=Avl(args[0], aircraft)
    avl.create_avl_file()
    avl.run_avl(args[1])

    
    
if __name__=="__main__":
    main(sys.argv[1], sys.argv[2])
