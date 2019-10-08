import io
import subprocess
import base64
from IPython.display import HTML

def make_movie(fname_glob, fname_out, fps_in=15, fps_out=15):
    # To force the frame rate of the input file (valid for raw formats only) to 1 fps and the frame rate of the output file to 24 fps:

    cmd = ['ffmpeg',
           '-y', # override existing file
#           '-r', str(fps_in),
           '-f', 'image2',
           '-framerate', str(fps_in),
           '-pattern_type', 'glob',
           '-i', fname_glob,
#           '-r', str(fps_out),
           '-pix_fmt', 'yuv420p',
           '-vcodec', 'libx264',
#           '-vf', 'scale=trunc(iw/2)*2:trunc(ih/2)*2',
           '-s','1200x800',
           '-f', 'mp4', fname_out]

    print('[make_mp4]: ffmpeg command:')
    print('{0:s}'.format(' '.join(cmd)))

    ret = subprocess.call(cmd)
    if ret == 0:
        print('[make_movie]: Successful execution. Output:')
        print('{0:s}'.format(fname_out))
    else:
        print('[make_movie]: subprocess.call returned {0:d}. Something went wrong.'.format(ret))
        
def display_movie(filename):

    video = io.open(filename, 'r+b').read()
    encoded = base64.b64encode(video)
    return HTML(data='''<video alt="test" controls>
                <source src="data:video/mp4;base64,{0}" type="video/mp4" />
                </video>'''.format(encoded.decode('ascii')))
