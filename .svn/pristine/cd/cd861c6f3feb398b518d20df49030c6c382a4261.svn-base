from pylab import xlim, errorbar, fill_between, array, append, plot, scatter, pcolor, show, figure, xlabel, ylabel, inf, sqrt, repeat, legend, xscale, yscale
import matplotlib.pyplot as plt
from numpy import cumsum, insert, zeros, sum
import numpy as np

import matplotlib
fig_xsize = 8
fig_ysize = 6

#matplotlib.rc('font', weight = 'bold')
matplotlib.rc('xtick', labelsize=14)#, labelweight = 'bold')
matplotlib.rc('ytick', labelsize=14)
matplotlib.rc('axes', labelsize = 14)
matplotlib.rc('ytick.minor', size=10)
matplotlib.rc('xtick.major', size=10)
matplotlib.rc("axes", linewidth=2.0, grid=True)#, labelweight = 'bold')
matplotlib.rc("lines", markeredgewidth=2.0, linewidth = 2)
#matplotlib.rc('figure', autolayout = True)
matplotlib.rc("grid", color='#888888', linewidth = 0.3)

from itertools import cycle
line_styles = ["-","--","-.",":"]

def apcolor(x= np.array([]), y=np.array([]), data=None, **kwargs):
    if x.size == 0:
        x = np.arange(0, data.shape[0]+1)
    if y.size == 0:
        y = np.arange(0, data.shape[1]+1)
    xpos = (x[1:] + x[:-1])/2.
    ypos = (y[1:] + y[:-1])/2.
    
    labels = []
    ax  = plt.pcolor(x, y, data, **kwargs)
    for yi in range(data.shape[0]):
        for xi in range(data.shape[1]):
            labels.append(plt.text(xpos[xi], ypos[yi], '%.2f' % data[yi, xi],
                                   horizontalalignment='center',
                                   verticalalignment='center',
                                   ))
    plt.xlim([x.min(), x.max()])
    plt.ylim([y.min(), y.max()])
    return ax, xpos, ypos

def errorMark(xaxis, bin_content,
              error = None,
              rel_error = None,
              color = 'blue',
              lw    = 2.5,
              label = None):

    if xaxis[-1] == inf:
        xbin = xaxis[-2] - xaxis[-3]
        xaxis[-1] = xaxis[-2] + xbin    

    xlim([xaxis[0], xaxis[-1]])

    xerror  = (xaxis[1:] - xaxis[:-1])/2
    error_xaxis = (xaxis[:-1]+xaxis[1:])/2
    
    if error == None:
        if rel_error == None:
            yerr=sqrt(bin_content)
        else:
            yerr = rel_error*bin_content
    else:
        yerr = error

    return errorbar(error_xaxis, bin_content, yerr=yerr, xerr = xerror, ls = 'none', ecolor = color, elinewidth = lw, capsize = 0, label = label)

def errorMarkVert(xaxis, bin_content,
                  yerror = None,
                  color = 'blue',
                  lw    = 2.5,
                  label = None):

    if xaxis[-1] == inf:
        xbin = xaxis[-2] - xaxis[-3]
        xaxis[-1] = xaxis[-2] + xbin    

    xlim([xaxis[0], xaxis[-1]])

    xerror  = (xaxis[1:] - xaxis[:-1])/2
    error_xaxis = (xaxis[:-1]+xaxis[1:])/2

    if yerror == None:
        yerr=sqrt(bin_content)
    else:
        yerr = yerror

    return errorbar(error_xaxis, bin_content, yerr=yerr, ls = 'none', ecolor = color, elinewidth = lw, capsize = lw*2, label = label)


def errorMarkHor(xaxis, bin_content,
                 lw = 2.5,
                 color = 'blue',
                 **kwargs):

    if xaxis[-1] == inf:
        xbin = xaxis[-2] - xaxis[-3]
        xaxis[-1] = xaxis[-2] + xbin    

    xlim([xaxis[0], xaxis[-1]])

    xerror  = (xaxis[1:] - xaxis[:-1])/2
    error_xaxis = (xaxis[:-1]+xaxis[1:])/2

    return errorbar(error_xaxis, bin_content, yerr=0, xerr = xerror, ls = 'none', elinewidth = lw, ecolor = color, capsize = 0, **kwargs)

def simpleMark(xaxis, bin_content,
               **kwargs):

    if xaxis[-1] == inf:
        xbin = xaxis[-2] - xaxis[-3]
        xaxis[-1] = xaxis[-2] + xbin    

    xlim([xaxis[0], xaxis[-1]])

    new_xaxis = (xaxis[:-1]+xaxis[1:])/2

    return plot(new_xaxis, bin_content, linestyle = 'None', **kwargs)


def weightedError(xaxis, bin_content,
                  bin_w2,
                  color = 'blue',
                  lw    = 2.5,
                  label = None):

    if xaxis[-1] == inf:
        xbin = xaxis[-2] - xaxis[-3]
        xaxis[-1] = xaxis[-2] + xbin    

    xlim([xaxis[0], xaxis[-1]])

    xerror  = (xaxis[1:] - xaxis[:-1])/2
    error_xaxis = (xaxis[:-1]+xaxis[1:])/2

    bin_error = sqrt(bin_w2)

    return errorbar(error_xaxis, bin_content, yerr=bin_error, xerr = None, ls = 'none', ecolor = color, elinewidth = lw, capsize = None, label = label)

def weightedError2(xaxis, bin_content,
                  bin_error,
                  color = 'blue',
                  lw    = 2.5,
                  label = None):

    if xaxis[-1] == inf:
        xbin = xaxis[-2] - xaxis[-3]
        xaxis[-1] = xaxis[-2] + xbin    

    xlim([xaxis[0], xaxis[-1]])

    xerror  = (xaxis[1:] - xaxis[:-1])/2
    error_xaxis = (xaxis[:-1]+xaxis[1:])/2


    return errorbar(error_xaxis, bin_content, yerr=bin_error, xerr = xerror, ls = 'none', ecolor = color, elinewidth = lw, capsize = 0, label = label)

def errorArea(xaxis, bin_up, bin_down,
              fillcolor = 'blue',
              hatch = None,
              hatchcolor = 'green',
              **kwargs):
    if hatch == None:
        hatchcolor = 'none'

    if xaxis[-1] == inf:
        xbin = xaxis[-2] - xaxis[-3]
        xaxis[-1] = xaxis[-2] + xbin

    xlim([xaxis[0], xaxis[-1]])

    fill_xaxis = repeat(xaxis,2)[1:-1]
    fill_bins_up  = repeat(bin_up,2)
    fill_bins_down = repeat(bin_down, 2)

    return fill_between(fill_xaxis , fill_bins_up, fill_bins_down, facecolor=fillcolor, edgecolor = hatchcolor, hatch = hatch,**kwargs)

def errorHatch(xaxis, bin_up, bin_down,
              hatch = '/',
              fillcolor = 'blue',
              **kwargs):

    if xaxis[-1] == inf:
        xbin = xaxis[-2] - xaxis[-3]
        xaxis[-1] = xaxis[-2] + xbin

    xlim([xaxis[0], xaxis[-1]])

    fill_xaxis = repeat(xaxis,2)[1:-1]
    fill_bins_up  = repeat(bin_up,2)
    fill_bins_down = repeat(bin_down, 2)

    return fill_between(fill_xaxis , fill_bins_up, fill_bins_down, edgecolor = hatchcolor, hatch = hatch, rasterized = True, **kwargs)


def unfilledBar(xaxis, bin_content, 
                errors = None,
                rel_errors = None,
                count_errors = False,
                color = 'blue',
                **kwargs):
    '''
    Produce a histogram bar plot a-la-ROOT. Returns plot axes.
    '''

    original_xaxis = xaxis
    if xaxis[-1] == inf:
        xbin = xaxis[-2] - xaxis[-3]
        xaxis[-1] = xaxis[-2] + xbin
    
    new_bins =  array(append(bin_content[0], bin_content))
    
    xlim([xaxis[0], xaxis[-1]])

    if rel_errors != None:
        errors = bin_content*rel_errors
    elif count_errors:
        errors = sqrt(bin_content)

    if errors != None:
        error_xaxis = array(xaxis)
        error_xaxis = error_xaxis[:-1] + (error_xaxis[1:] - error_xaxis[:-1])/2.
        errorbar(error_xaxis, bin_content, yerr = errors, color = color, fmt='.') 

    return plot(xaxis , new_bins, drawstyle = 'steps', color = color,  **kwargs)

def stackedBarPlot(xaxis, hist_list,
                   colors = ['r', 'b', 'g', 'm', 'c', 'y'],
                   linewidth = 1,
                   fill = True,
                   outer_line = True,
                   hatch = [None]*6,
                   zorder = -9,
                   data_axis = 0,
                   relative = False,
                   **kwargs):

    bin_array = array(hist_list)
    cumsum_axis = 0
    # Do the cumulative sum of the data
    bins_cumsum = insert(cumsum(bin_array, axis = cumsum_axis), 0, zeros(len(xaxis)-1), cumsum_axis)
    if relative:
        bins_cumsum /= bins_cumsum[-1, :]

    # Fix the xaxis
    original_xaxis = xaxis
    if xaxis[-1] == inf:
        xbin = xaxis[-2] - xaxis[-3]
        xaxis[-1] = xaxis[-2] + xbin

    xlim([xaxis[0], xaxis[-1]])
    fill_xaxis = repeat(xaxis,2)[1:-1]


    areas = []
    lines = []
    rectangles = []
    print 'Cumsum results ', sum(bins_cumsum, axis = 1)

    nhistograms = bins_cumsum.shape[cumsum_axis]

    for index in range(1, bins_cumsum.shape[cumsum_axis]):
        if outer_line:
            lines.append(unfilledBar(original_xaxis, bins_cumsum[index, :], color = colors[index-1], linewidth = linewidth))
            lwbox = 1
        else:
            lwbox = 0
            #lines.append(unfilledBar(original_xaxis, bins_cumsum[:, index], color = colors[index-1], label = labels[index-1]))

        if hatch[index-1]:
            facecolor = 'white'
        else:
            facecolor = colors[index-1]
        if fill:
            areas.append( fill_between(repeat(xaxis,2)[1:-1], repeat(bins_cumsum[index,:],2), repeat(bins_cumsum[index-1,:],2),
                                       zorder = zorder-nhistograms+index, linewidth = 3,
                                       facecolor=facecolor, edgecolor = colors[index-1], hatch = hatch[index-1], **kwargs ))

            rectangles.append(matplotlib.patches.Rectangle((0,0), 1, 1, fc = facecolor, linewidth=2, zorder = zorder-1,
                                                           edgecolor = colors[index-1], hatch = hatch[index-1]))


    return lines, areas, rectangles
    # Use legend(rectangles, labels)

# def sampleContribution(xaxis, hist_list,
#                        colors = ['r', 'b', 'g', 'm', 'c', 'y'],
#                        linewidth = 1,
#                        fill = True,
#                        outer_line = True,
#                        hatch = [None]*6,
#                        zorder = -9,
#                        data_axis = 0,
#                        **kwargs):

#     bin_array = array(hist_list)
#     cumsum_axis = 0
#     # Do the cumulative sum of the data
#     bins_cumsum = insert(cumsum(bin_array, axis = cumsum_axis), 0, zeros(len(xaxis)-1), cumsum_axis)
#     bins_cumsum /= 

#     # Fix the xaxis
#     original_xaxis = xaxis
#     if xaxis[-1] == inf:
#         xbin = xaxis[-2] - xaxis[-3]
#         xaxis[-1] = xaxis[-2] + xbin

#     xlim([xaxis[0], xaxis[-1]])
#     fill_xaxis = repeat(xaxis,2)[1:-1]


#     areas = []
#     lines = []
#     rectangles = []
#     print 'Cumsum results ', sum(bins_cumsum, axis = 1)

#     nhistograms = bins_cumsum.shape[cumsum_axis]

#     for index in range(1, bins_cumsum.shape[cumsum_axis]):
#         if outer_line:
#             lines.append(unfilledBar(original_xaxis, bins_cumsum[index, :], color = colors[index-1], linewidth = linewidth))
#             lwbox = 1
#         else:
#             lwbox = 0
#             #lines.append(unfilledBar(original_xaxis, bins_cumsum[:, index], color = colors[index-1], label = labels[index-1]))

#         if hatch[index-1]:
#             facecolor = 'white'
#         else:
#             facecolor = colors[index-1]
#         if fill:
#             areas.append( fill_between(repeat(xaxis,2)[1:-1], repeat(bins_cumsum[index,:],2), repeat(bins_cumsum[index-1,:],2),
#                                        zorder = zorder-nhistograms+index, linewidth = lwbox,
#                                        facecolor=facecolor, edgecolor = colors[index-1], hatch = hatch[index-1], **kwargs ))

#             rectangles.append(matplotlib.patches.Rectangle((0,0), 1, 1, fc = facecolor, linewidth=lwbox, zorder = zorder-1,
#                                                            edgecolor = colors[index-1], hatch = hatch[index-1]))


#     return lines, areas, rectangles
