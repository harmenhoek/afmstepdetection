# AFM step detector

# Imports
import numpy as np
from matplotlib import cm
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy import interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import optimize
import json
import csv
import os
from datetime import datetime


fntsize = 20
font = {'size': fntsize}
plt.rc('font', **font)
mpl.rcParams['axes.linewidth'] = 2 #set the value globally


# Settings
width = 512
height = 256

filename = r'data/Example1.spm'

filename_output = 'output'

r_start = 0.5
r_end = 1.5
interpol_stepsize = 1

custom_theta = False
theta_set = [0, 10, 20]
theta_steps = 20

# horizontal fit settings
fit_offset = 5  # in resolution aka datapoints
residuals_threshold = 20000

# arrow for theta indication
R_arrow_out = 10  # how far the arrow starts beyond r_end




# Initialize settings
if custom_theta:
    theta_range = theta_set
else:
    theta_range = np.arange(0, 360, 360/theta_steps)

# Functions
def data_import(filename):
    data_start = 'end of header'
    data_end = 'end of experiment'
    with open(filename, 'r') as f:
        data = f.read()
    # remove header and footer of data
    data = data.split(data_start, 1)[1]
    data = data.split(data_end, 1)[0]
    # convert data str to list
    data = list(data.split("\n"))
    # remove empty values at start and end
    while data[0] == '':
        del data[0]
    while data[-1] == '':
        del data[-1]
    data = [float(i) for i in data]  # convert to floats
    # verify size
    if len(data) != width * height:
        raise Exception("Number of datapoints is not corresponding to width*height.")

    return data



def calc_R(x,y, xc, yc):
    """ calculate the distance of each 2D points from the center (xc, yc) """
    return np.sqrt((x-xc)**2 + (y-yc)**2)

def f(c, x, y):
    """ calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc) """
    Ri = calc_R(x, y, *c)
    return Ri - Ri.mean()

def leastsq_circle(x, y):
    # coordinates of the barycenter
    x_m = np.mean(x)
    y_m = np.mean(y)
    center_estimate = x_m, y_m
    center, ier = optimize.leastsq(f, center_estimate, args=(x, y))
    xc, yc = center
    Ri = calc_R(x, y, *center)
    R = Ri.mean()
    residu = np.sum((Ri - R)**2)
    return xc, yc, R, residu



def tellme(s):
    print(s)
    plt.title(s, fontsize=16)
    plt.draw()

def get_circle_v2():
    plt_width = 20
    plt_interpolation = 'nearest'  # 'off'

    fig, axes = plt.subplots(1, 1, figsize=(plt_width, height / width * plt_width))
    cs = plt.imshow(np.flipud(data), interpolation=plt_interpolation, cmap=plt.get_cmap('YlOrBr_r'),
                    vmin=np.min(data), vmax=np.max(data), origin='lower')
    plt.colorbar(cs)

    pts, x, y = [], [], []
    while True:

        if len(pts) < 3:
            while len(pts) < 3:
                tellme('Select 3 points with mouse')
                pts = np.asarray(plt.ginput(3, timeout=-1))
                if len(pts) < 3:
                    tellme('Too few points, starting over')
            x = pts[:, 0]
            y = pts[:, 1]

        xc, yc, R, residu = leastsq_circle(x, y)
        circle = plt.Circle((xc, yc), R, fill=False, linewidth=5)
        axes.add_artist(circle)
        plt.draw()
        p, = plt.plot(x, y, 'o', color='black')
        tellme(f'Current residual: {round(residu, 2)}. Mouse click to select an extra point, key press to fit.')
        if plt.waitforbuttonpress():
            break
        pt = np.asarray(plt.ginput(n=1, timeout=-1))
        x = np.append(x, pt[0, 0])
        y = np.append(y, pt[0, 1])
        p.remove()
        circle.remove()
    return xc, yc, R

def get_circle_v1():
    plt_width = 20
    plt_interpolation = 'nearest'  # 'off'

    fig, axes = plt.subplots(1, 1, figsize=(plt_width, height / width * plt_width))
    cs = plt.imshow(np.flipud(data), interpolation=plt_interpolation, cmap=plt.get_cmap('YlOrBr_r'),
                    vmin=np.min(data), vmax=np.max(data), origin='lower')
    plt.colorbar(cs)
    while True:
        pts = np.asarray(plt.ginput(n=2, timeout=-1))
        draw_circle = plt.Circle((pts[0]), np.linalg.norm(pts[0]-pts[1]), fill=False)
        axes.add_artist(draw_circle)
        p, = plt.plot([pts[0][0]], [pts[0][1]], 'o')
        tellme('Happy? Key click for yes, mouse click for no')
        if plt.waitforbuttonpress():
            break
        p.remove()
        draw_circle.remove()
    plt.close(fig)
    return pts[0][0], pts[0][1], np.linalg.norm(pts[0]-pts[1])

def convert_coordinates(x0_c, y0_c, x, y):
    # x0_c, y0_c are center coordinates of circle in x0, y0 system
    # x, y are the coordinates to convert from x0,y0 system to x0_c, y0_c system
    x_c = x + x0_c
    y_c = y + y0_c
    return x_c, y_c

def convert_coordinates_back(x0_c, y0_c, x_c, y_c):
    # x0_c, y0_c are center coordinates of circle in x0, y0 system
    # x, y are the coordinates to convert from x0,y0 system to x0_c, y0_c system
    x = x_c - x0_c
    y = y_c - y0_c
    return x, y

def cart2pol(x, y, output='rad'):
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    if output == 'deg':
        theta = np.rad2deg(theta)
    return(r, theta)

def pol2cart(r, theta, input='rad'):
    if input == 'deg':
        theta = np.deg2rad(theta)
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return(x, y)

def verify_coordinates(x, y):
    # check if coordinates are inside image, otherwise ignore
    check = []
    for xi, yi in zip(x, y):
        if xi < 0 or xi > width or yi < 0 or yi > height:
            check.append(False)
        else:
            check.append(True)
    return check

def fithor(x, y, x_threshold, part):
    if part == 'below':
        r_fit = [xi for xi in x if xi < (x_threshold)]
        z_fit = [yi for xi, yi in zip(x, y) if xi < (x_threshold)]
    elif part == 'above':
        r_fit = [xi for xi in x if xi > (x_threshold)]
        z_fit = [yi for xi, yi in zip(x, y) if xi > (x_threshold)]
    else:
        raise Exception('No valid part selected. Try below or above.')
    [coefficients, residuals, _, _, _] = np.polyfit(r_fit, z_fit, 0, full=True)
    return r_fit, np.poly1d(coefficients), residuals

def namestr(obj, namespace):
    ''' To print the var name inside def save_csv '''
    return [name for name in namespace if namespace[name] is obj]

def save_csv(filename, data, header=[], headername='radii'):
    '''
    This def saves data in form of dict to filename. Basically 2 types:
    1 --> where each key in dict has 1 value or 1 list (not np!). key will be on column 1, data following in the columns after
        example:
            output = {
                'number_of_experiments': a1,
                'radii': [r1, r2, r3, r4],
                'bananas': [b1, b2],
            }
        output save_csv(filename, output) gives:
            number_of_experiments, a1
            radii, r1, r2, r3, r4
            bananas, b1, b2
    2 --> one more nested list. so dict-->key-->list-->list. In this key each main list will be saved on different rows
        the start and end of data is indicated with string. Column 1 will be header.
        example:
            raw_data = {
                'data': [
                            [a1,a2,a3],
                            [b1,b2,b3],
                ],
            }
        output save_csv(filename, raw_data, header=radii, headername='Radii'), where radii=[r1,r2] gives:
            START OF raw_data DATA.
            COLUMN 1 SHOWS radii (Radii)
            r1, a1, a2, a3
            r2, b1, b2, b3
            END OF raw_data DATA
    Combination of the 2 is also possible.
    '''
    with open(filename, 'w') as f:
        f.write(f'START OF HEADER\n')
        f.write(
            f'This output file was auto generated using {os.path.basename(__file__)}\nDate time: {datetime.now().strftime("%d/%m/%Y %H:%M:%S")}\n')
        f.write(f'SETTINGS\n')
        f.write(f'width,{width}\nheight,{height}\nfilename,{filename}\nr_start,{r_start}\nr_end,{r_end}\ninterpol_stepsize,{interpol_stepsize}\nfit_offset,{fit_offset}\nresiduals_threshold,{residuals_threshold}\n')
        f.write(f'END OF HEADER\nSTART OF DATA\n')

        for key in data.keys():
            if type(data[key]) == list:  # type 2
                if type(data[key][0]) == list: # raw data
                    f.write(f'START OF {key} DATA. COLUMN 1 SHOWS r_range (rawdata_x)\n ')# \nCOLUMN 1 SHOWS {namestr(header, globals())[0]} ({headername})\n') TODO temp fix for Tkinter
                    for item, hdr in zip(data[key], header):
                        f.write(f'{hdr},')
                        for subitem in item:
                            f.write(f'{subitem},')
                        f.write(f'\n')
                    f.write(f'END OF {key} DATA \n')
                else:  # type 1
                    f.write(f'{key},')
                    for item in data[key]:
                        f.write(f'{item},')
                    f.write(f'\n')
            else:
                f.write(f'{key},{data[key]} \n')
        f.write(f'END OF DATA')


import tkinter
import tkinter as tk
from tkinter import Frame, Tk, Button

from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
# Implement the default Matplotlib key bindings.
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure

import numpy as np


root = tkinter.Tk()
root.title("Embedding in Tk")
# root.iconbitmap('favicon.ico')

# # general canvas with grid
# cnv = tk.Canvas(root, width=800, height=400, background='white')
# # cnv.grid(columnspan=2, rowspan=2)
#
#
# # generate figure
# fig = Figure(figsize=(5, 4), dpi=100)
# # t = np.arange(0, 3, .01)
# # fig.add_subplot(111).plot(t, 2 * np.sin(2 * np.pi * t))
#
# canvas = FigureCanvasTkAgg(fig, master=root)  # A tk.DrawingArea.
# canvas.draw()
# canvas.get_tk_widget().grid(column=0, row=0)
# # canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
#
# # toolbar for figure
# frame = Frame(root)
# frame.grid(row=2, column=0, sticky='we')
# toolbar = NavigationToolbar2Tk(canvas, frame)
# toolbar.update()
# # canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
# canvas.get_tk_widget().grid(column=0, row=1)

# def on_key_press(event):
#     print("you pressed {}".format(event.key))
#     key_press_handler(event, canvas, toolbar)
#
# canvas.mpl_connect("key_press_event", on_key_press)


def _quit():
    root.quit()     # stops mainloop
    root.destroy()  # this is necessary on Windows to prevent
                    # Fatal Python Error: PyEval_RestoreThread: NULL tstate

# quit button
# button = tkinter.Button(master=root, text="Quit", command=_quit)
# button.pack(side=tkinter.BOTTOM)
# button.grid(column=0, row=2)

## ----------------- RUN CODE -----------------

def run_afmstepdetector():

    filename_outputcsv = filename_output + '.csv'
    filename_outputjson = filename_output + '.json'

    global data
    data = data_import(filename)
    data = np.reshape(data, (height, width)) * 1e9

    # Get GUI to select the center of the droplet and the outside border
    # x0_c, y0_c, R = get_circle_v1()
    x0_c, y0_c, R = get_circle_v2()

    # R = 85.79
    # x0_c, y0_c = 260, 110

    print(x0_c, y0_c, R)

    # x0, y0 are at the left-bottom corner of the image.
    # x0_c, y0_c are at the center of the drop

    # Approach
    # 1. Define the line for interpolation (datapoints = stepsize of image)
    # 2. Interpolate the surface
    # 3. Evaluate the interpolated surface at the data points of the line

    # radii to evaluate (0 = center of drop)
    r_range = np.arange(r_start * R, r_end * R, interpol_stepsize)

    # get colormap for figures
    clrmp1 = cm.get_cmap('jet', 12)
    clrmp = clrmp1(np.linspace(0, 1, theta_steps))

    # interpolate surface (resolution x and y is 1)
    f = interpolate.interp2d(np.arange(0, width), np.flipud(np.arange(0, height)), data)

    # get basic plots ready
    plt_width = 20
    plt_interpolation = 'nearest'  # 'off'
    fig, axes = plt.subplots(1, 1, figsize=(plt_width, height / width * plt_width))
    fig2, axes2 = plt.subplots(1, 1, figsize=(plt_width, height / width * plt_width))
    cs = axes2.imshow(np.flipud(data), interpolation=plt_interpolation, cmap=plt.get_cmap('YlOrBr_r'),
                      vmin=np.min(data), vmax=np.max(data), origin='lower')
    divider = make_axes_locatable(axes2)
    cax = divider.append_axes("right", size="2%", pad=0.05)
    cbar = plt.colorbar(cs, cax=cax)
    cbar.ax.get_yaxis().labelpad = 15
    cbar.ax.set_ylabel('Height [nm]', rotation=90)
    draw_circle = plt.Circle((x0_c, y0_c), R, fill=False, linewidth=4)
    axes2.add_artist(draw_circle)
    axes2.plot([x0_c], [y0_c], 'o', markersize=15, color='black')

    output = {
        'step_heights': [],
        'fit_low_residuals': [],
        'fit_upp_residuals': [],
        'fit_low_nrdatapoints': [],
        'fit_upp_nrdatapoints': [],
        'ignored': [],
        'stepheight_median': [],
        'stepheight_average': [],
        'stepheight_std': [],
        'stepheight_var': [],
        'angles': theta_range.tolist(),
        'number_angles_total': len(theta_range),
        'number_angles_afteromit': [],
        'rawdata_x': r_range.tolist(),
        'rawdata_y': [],
        'plateau_xc_yx_R': [x0_c, y0_c, R],
    }

    # iterate over theta around the step
    for idx, theta in enumerate(theta_range):
        r = r_range
        # get the cartesian coordinates corresponding to current r(s) and theta
        x, y = pol2cart(r, theta, input='deg')
        x_c, y_c = convert_coordinates(x0_c, y0_c, x, y)
        # verify if coordinates are not outside image, if so remove those
        check = verify_coordinates(x_c, y_c)
        x_c = np.array(x_c)[check].tolist()
        y_c = np.array(y_c)[check].tolist()
        x = np.array(x)[check].tolist()
        y = np.array(y)[check].tolist()
        r = r[check].tolist()

        # evaluate interpolated surface the coordinates
        z = np.array([f(xi, yi)[0] for xi, yi in zip(x_c, y_c)])

        # fit horizontal lines tp top and bottom part
        r_part1, fit1, residuals1 = fithor(r, z, R + fit_offset, 'above')  # lower
        r_part2, fit2, residuals2 = fithor(r, z, R - fit_offset, 'below')  # upper
        z_offset = fit1(r_part1)[0]

        output['rawdata_y'].append(z.tolist())
        output['step_heights'].append(fit2(r_part2)[0] - fit1(r_part1)[0])
        output['fit_low_residuals'].append(residuals1.tolist()[0])
        output['fit_upp_residuals'].append(residuals2.tolist()[0])
        output['fit_low_nrdatapoints'].append(len(r_part1))
        output['fit_upp_nrdatapoints'].append(len(r_part2))

        # only if fit is decent, plots can be shown
        if residuals1 < residuals_threshold and residuals2 < residuals_threshold:
            output['ignored'].append(False)
            axes.plot(r_part1, fit1(r_part1) - z_offset, ':', color=clrmp[idx], linewidth=5)
            axes.plot(r_part2, fit2(r_part2) - z_offset, ':', color=clrmp[idx], linewidth=5)
            axes.plot(r, z - z_offset, '.-', markersize=10, linewidth=3, color=clrmp[idx])
            axes2.plot(x_c, y_c, '.-', markersize=2, color=clrmp[idx])

            x_c = [xi for xi, ri in zip(x_c, r_range) if ri < (R - fit_offset) or ri > (R + fit_offset)]
            y_c = [yi for yi, ri in zip(y_c, r_range) if ri < (R - fit_offset) or ri > (R + fit_offset)]
            axes2.plot(x_c, y_c, '.-', markersize=12, color=clrmp[idx])
        else:
            output['ignored'].append(True)
            print(
                f'Angle theta={theta} was omitted, since one of the residuals exceeded the set threshold ({round(residuals1[0])}, {round(residuals2[0])} > {residuals_threshold}).')

    output['stepheight_average'] = np.mean(output['step_heights'])
    output['stepheight_median'] = np.median(output['step_heights'])
    output['stepheight_std'] = np.std(output['step_heights'])
    output['stepheight_var'] = np.var(output['step_heights'])
    output['number_angles_afteromit'] = np.size(output['ignored']) - np.count_nonzero(output['ignored'])

    save_csv(filename_outputcsv, output, header=r_range, headername='rawdata_x')
    output_json = json.dumps(output, sort_keys=True, indent=4)

    # Plot arrow for theta angle indication
    style = "Simple, tail_width=4, head_width=15, head_length=12"
    kw = dict(arrowstyle=style, color="k", hatch='|')
    (x, y) = pol2cart(r_end * R + R_arrow_out, 25, input='deg')
    (x, y) = convert_coordinates(x0_c, y0_c, x, y)
    a3 = patches.FancyArrowPatch((x0_c + r_end * R + R_arrow_out, y0_c), (x, y), connectionstyle="arc3,rad=.15", **kw)
    axes2.add_patch(a3)
    axes2.plot([x0_c, x0_c + r_end * R + R_arrow_out + 10], [y0_c, y0_c], ':', color='black', zorder=1, linewidth=5)
    # print text theta
    (x, y) = pol2cart(r_end * R + R_arrow_out + 5, 25 / 2, input='deg')
    (x, y) = convert_coordinates(x0_c, y0_c, x, y)
    axes2.annotate(r'$\theta$', (x, y), fontsize=fntsize, va='center')

    axes.set_xlabel('Radial distance from center circle [pixels]')
    axes.set_ylabel('Height [nm]')
    axes2.set_xlabel('y [pixels]')
    axes2.set_ylabel('x [pixels]')

    axes.set_xlim([r_range[0], r_range[-1:]])
    axes.plot((r_range[0], r_range[-1]), (0, 0), '-', color='black', zorder=1, linewidth=3)

    sm = plt.cm.ScalarMappable(cmap=cm.get_cmap('jet', 360), norm=plt.Normalize(vmin=min(theta_range), vmax=360))
    divider = make_axes_locatable(axes)
    cax = divider.append_axes("right", size="2%", pad=0.05)
    cbar = plt.colorbar(sm, ticks=range(0, 361, 45), spacing='proportional', cax=cax)  # , boundaries=theta_range)
    cbar.ax.get_yaxis().labelpad = 15
    cbar.ax.set_ylabel(r'$\theta$ [deg]', rotation=90)

    plt.show()




def run():
    # get all input buttons
    global filename, filename_output, width, height, r_start, r_end, theta_steps, custom_theta, theta_set, fit_offset, residuals_threshold
    filename = e_filename.get()
    filename_output = e_filenameout.get()
    width = int(e_width.get())
    height = int(e_height.get())
    r_start = float(e_rstart.get())
    r_end = float(e_rend.get())
    theta_steps = int(e_thetasteps.get())
    custom_theta = bool(e_customtheta_value.get())
    theta_set = list(e_thetaset.get().split())
    theta_set = [float(i) for i in theta_set]
    fit_offset = int(e_fitoffset.get())
    residuals_threshold = float(e_resthreshold.get())

    print('filename: ', filename, f"({type(filename)})")
    print('filename_output: ', filename_output, f"({type(filename_output)})")
    print('width: ', width, f"({type(width)})")
    print('height: ', height, f"({type(height)})")
    print('r_start: ', r_start, f"({type(r_start)})")
    print('r_end: ', r_end, f"({type(r_end)})")
    print('theta_steps: ', theta_steps, f"({type(theta_steps)})")
    print('custom_theta: ', custom_theta, f"({type(custom_theta)})")
    print('theta_set: ', theta_set, f"({type(theta_set)})")
    print('fit_offset: ', fit_offset, f"({type(fit_offset)})")
    print('residuals_threshold: ', residuals_threshold, f"({type(residuals_threshold)})")
    # print('r_start: ', r_start, f"({type(r_start)})")
    # exit()

    # run code as usual
    run_afmstepdetector()
    return


## ----------------- RIGHT SIDE OF WINDOW -----------------

frame_right = Frame(root, width=400)
frame_right.grid(row=1, column=1, columnspan=1, sticky='nwe')

r = 0
tk.Label(frame_right, text="AFM Step Detector V1", font=('Arial', 14)).grid(row=r, column=0, columnspan=3, pady=10)
r += 1
tk.Label(frame_right, text='Settings', font=('Arial', 12)).grid(row=r, column=0, padx=2, pady=5, sticky='w')

r += 1
l_filename = tk.Label(frame_right, text='Input filename').grid(row=r, column=0, sticky='e')
e_filename = tk.Entry(frame_right)
e_filename.insert(0, filename)
e_filename.grid(row=r, column=1)

r += 1
l_filenameout = tk.Label(frame_right, text='Output filename').grid(row=r, column=0, sticky='e')
e_filenameout = tk.Entry(frame_right)
e_filenameout.insert(0, filename_output)
e_filenameout.grid(row=r, column=1)
l_filenameout2 = tk.Label(frame_right, text='.csv / .json').grid(row=r, column=2, sticky='w')

r += 1
l_width = tk.Label(frame_right, text='Image width').grid(row=r, column=0, sticky='e')
e_width = tk.Entry(frame_right)
e_width.insert(0, width)
e_width.grid(row=r, column=1)
l_width2 = tk.Label(frame_right, text='pixels').grid(row=r, column=2, sticky='w')

r += 1
l_height = tk.Label(frame_right, text='Image height').grid(row=r, column=0, sticky='e')
e_height = tk.Entry(frame_right)
e_height.insert(0, height)
e_height.grid(row=r, column=1)
l_height2 = tk.Label(frame_right, text='pixels').grid(row=r, column=2, sticky='w')

r += 1
l_rstart = tk.Label(frame_right, text='r_start').grid(row=r, column=0, sticky='e')
e_rstart = tk.Entry(frame_right)
e_rstart.insert(0, r_start)
e_rstart.grid(row=r, column=1)
l_rstart2 = tk.Label(frame_right, text='R').grid(row=r, column=2, sticky='w')

r += 1
l_rend = tk.Label(frame_right, text='r_end').grid(row=r, column=0, sticky='e')
e_rend = tk.Entry(frame_right)
e_rend.insert(0, r_end)
e_rend.grid(row=r, column=1)
l_rend2 = tk.Label(frame_right, text='R').grid(row=r, column=2, sticky='w')

r += 1
l_thetasteps = tk.Label(frame_right, text='Theta step size').grid(row=r, column=0, sticky='e')
e_thetasteps = tk.Entry(frame_right)
e_thetasteps.insert(0, theta_steps)
e_thetasteps.grid(row=r, column=1)
l_thetasteps2 = tk.Label(frame_right, text='degrees').grid(row=r, column=2, sticky='w')

r += 1
e_customtheta_value = tk.BooleanVar()
e_customtheta_value.set(False)
e_customtheta = tk.Checkbutton(frame_right, text='Custom theta', variable=e_customtheta_value).grid(row=r, column=1, sticky='w')

r += 1
l_thetaset = tk.Label(frame_right, text='Custom theta range').grid(row=r, column=0, sticky='e')
e_thetaset = tk.Entry(frame_right)
e_thetaset.insert(0, theta_set)
e_thetaset.grid(row=r, column=1)
l_thetaset2 = tk.Label(frame_right, text='degrees').grid(row=r, column=2, sticky='w')

r += 1
l_fitoffset = tk.Label(frame_right, text='Fit offset').grid(row=r, column=0, sticky='e')
e_fitoffset = tk.Entry(frame_right)
e_fitoffset.insert(0, fit_offset)
e_fitoffset.grid(row=r, column=1)
l_fitoffset2 = tk.Label(frame_right, text='resolution (pixels)').grid(row=r, column=2, sticky='w')

r += 1
l_resthreshold = tk.Label(frame_right, text='Residual threshold').grid(row=r, column=0, sticky='e')
e_resthreshold = tk.Entry(frame_right)
e_resthreshold.insert(0, residuals_threshold)
e_resthreshold.grid(row=r, column=1)
l_resthreshold2 = tk.Label(frame_right, text='(above discarded)').grid(row=r, column=2, sticky='w')

# r += 1
# tk.Label(frame_right, text='Run', font=('Arial', 12)).grid(row=r, column=0, padx=2, pady=5, sticky='w')

r += 1
b_run = Button(frame_right, text='RUN', padx=20, pady=10, command=run).grid(row=r, column=0, columnspan=3, pady=10)
r += 1
b_quit = tkinter.Button(master=root, text="Quit", command=_quit).grid(row=r, column=0, columnspan=3)






tkinter.mainloop()
# If you put root.destroy() here, it will cause an error if the window is
# closed with the window manager.