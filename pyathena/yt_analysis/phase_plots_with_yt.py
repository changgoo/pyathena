import yt
import pyathena.yt_analysis.ytathena as ya

tigress_unit_system=yt.UnitSystem('tigress','pc','Msun','Myr',)
tigress_unit_system['velocity']='km/s'
tigress_unit_system['magnetic_field']='uG'

bin_fields=[]

bin_fields.append(['nH','pok'])
bin_fields.append(['nH','temperature'])
bin_fields.append(['velocity_magnitude','sound_speed'])
bin_fields.append(['nH','velocity_magnitude'])
bin_fields.append(['radius','velocity_magnitude'])
bin_fields.append(['radius','temperature'])
bin_fields.append(['radius','nH'])
bin_fields.append(['nH','magnetic_field_magnitude'])
bin_fields.append(['temperature','magnetic_field_magnitude'])
bin_fields.append(['velocity_magnitude','magnetic_field_magnitude'])
bin_fields.append(['kinetic_energy','magnetic_energy'])

extrema={}
extrema['nH']=(1.e-4,1.e4)
extrema['pok']=(1,1.e8)
extrema['temperature']=(1,1.e9)
extrema['velocity_magnitude']=(1.e-2,1.e4)
extrema['sound_speed']=(1.e-2,1.e4)
extrema['radius']=(0,32)
extrema['mag_pok']=(1,1.e8)
extrema['magnetic_field_magnitude']=(1.e-2,1.e4)
extrema['kinetic_energy']=(1.e38,1.e48)
extrema['magnetic_energy']=(1.e38,1.e48)

logs={}
logs['radius']=False

field='cell_mass'

def phase(filename,out_dir='',write_file=True,write_figure=True):
    ds=yt.load(filename,units_override=ya.unit_base,unit_system=tigress_unit_system)
    ya.add_yt_fields(ds,cooling=True,mhd=True,rotation=False)

    sp=ds.sphere(ds.domain_center,ds.domain_right_edge[0])

    total_mass=sp.sum('cell_mass')
    print total_mass.in_units('Msun')

    for bf in bin_fields:
        xbin,ybin=bf
        pdf=yt.create_profile(sp,bf,field,extrema=extrema,logs=logs,
                              n_bins=128,weight_field=None)
        outhead1='{}{}_{}_{}_{}'.format(out_dir,ds,xbin,ybin,field)
        if write_file: pdf.save_as_dataset(outhead1)
        if write_figure: 
            p=pdf.plot()
            outhead2='{}.png'.format(outhead1)
            p.save(outhead2)

def phase_plot(phase_dataset,out_dir=''):
    pdf_ds=yt.load(phase_dataset)
    data=pdf_ds.data
    p=yt.PhasePlot(data,pdf_ds.x_field,pdf_ds.y_field,'cell_mass',weight_field=None)
    p.save('{}{}'.format(out_dir,pdf_ds).replace('.h5','.png'))
