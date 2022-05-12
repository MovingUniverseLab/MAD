import os
import pandas as pd
import numpy as np
from flask import Flask, Response, flash, request, redirect, url_for, send_from_directory, render_template
from werkzeug.utils import secure_filename
import sqlite3
from sqlalchemy import create_engine, text
from sqlalchemy.sql import select
import datetime
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import io
import base64
import numexpr as ne
from flask import make_response


app = Flask(__name__)


# The database.
engine = create_engine('sqlite:///microlensing.db')


@app.route('/', methods=['GET', 'POST'])
def query_db():
    """
    Homepage. Provides the interface to query the database.
    """
    # df_query_result is the result of the SQL query performed. 
    # It is passed around the whole app (global variable),
    # and gets overwritten each time there is a new query.
    global df_query_result

    if request.method == 'POST':
        # Use pandas to perform the SQL query, then convert it to html
        # so we can pass the results into a table for display on the page.
        query_str = request.form['query']
        with engine.connect() as conn:
            df_query_result = pd.read_sql(query_str, conn)
            df_html = df_query_result.to_html()
            
            # If 'alert_name' is returned in the query, provide the
            # option to view the lightcurves in the template.
            # Otherwise it doesn't get shown.
            if 'alert_name' in df_query_result.columns:
                browse_lc=True
            else:
                browse_lc=False

        # Display the results (if there are any).
        if len(df_query_result) == 0:
            return render_template('display_empty.html', 
                                   query_str=query_str,
                                   query_db=url_for('query_db'))
        else:        
            return render_template('display_table.html', 
                                   query_str=query_str,
                                    html_table=df_html,        
                                   download_csv=url_for('download_csv', query_str=query_str),
                                   browse_lc=browse_lc,
                                   browse_lightcurves=url_for('browse_lightcurves'))
        
    return render_template('query.html')


@app.route('/download_csv/<query_str>', methods=['GET', 'POST'])
def download_csv(query_str):
    """
    Download the result of the SQL query as a CSV file.
    Note: can read the result with pandas read_csv and it will
    give you the column names.
    """
    with engine.connect() as conn:
        df = pd.read_sql(query_str, conn)
    resp = make_response(df.to_csv())
    resp.headers["Content-Disposition"] = "attachment; filename=export.csv"
    resp.headers["Content-Type"] = "text/csv"
    
    return resp


@app.route('/browse_lightcurves', methods=['GET', 'POST'])
def browse_lightcurves():
    """
    Page that lists all the lightcurve from the query.
    Each entry links to a page that shows the lightcurve.
    """
    # Get a unique list of the alert names.
    items = df_query_result['alert_name'].tolist()
    alert_names = list(dict.fromkeys(items))
    
    return render_template('lightcurves_list.html', 
                            alert_names=alert_names)


@app.route('/fig/<alert_name>')
def plot_lightcurve(alert_name):
    """
    Page that shows the MOA lightcurve in magnitude space.
    """
    # Grab the hjd, mag, mag_err corresponding to the alert name.
    query_str = 'SELECT hjd, mag, mag_err FROM photometry WHERE alert_name = "' + alert_name + '"' 
    db_info = engine.execute(query_str).fetchall()
    time = np.array([info[0] for info in db_info])
    mag = np.array([info[1] for info in db_info])
    mag_err = np.array([info[2] for info in db_info])
    
    #####
    # We don't want to plot the whole time series (can be 10+ years).
    # So we are only going to grab the data from the year of the alert
    # and the year preceeding it.
    #####
    # Get the year of the alert.
    YY = alert_name[2:4]
    
    # The *1 is just a dumb trick to turn it into an integer.
    year = ne.evaluate(YY) * 1 
    
    # HJD date corresponding to 1 January 2000.
    hjd_jan_00 = 1154

    # Calculate the dates that correspond to a time span of the 
    # year of the alert and the year before it.
    # Why is it year and year+2? Should be year -1 and year+1?
    start_date = hjd_jan_00 + 365.25 * (year)
    end_date = hjd_jan_00 + 365.25 * (year + 2) # This gives us data through that year's alert season.
    
    # Now only keep things from the year of and before.
    keep_idx = np.where((time < end_date) & (time > start_date))[0]
    
    # Finally... make the plot.
    fig = create_figure(time[keep_idx], mag[keep_idx], mag_err[keep_idx], alert_name)
    
    # Get a unique list of the alert names.
    items = df_query_result['alert_name'].tolist()
    alert_names = list(dict.fromkeys(items))
    
    # Figure out which entry in the list this is so we know which template to use below.
    # Tried to get all those if statement below into the template but couldn't get it quite to work...
    n_lc = len(alert_names)
    ii = alert_names.index(alert_name)

    # Catch edge case with only one result.
    if n_lc == 1:
        return render_template('show_lightcurve_one.html', 
                                home=url_for('query_db'),
                                image=fig)
    # First page.
    if ii == 0:
        return render_template('show_lightcurve_first.html', 
                                home=url_for('query_db'),
                                next_page=url_for('plot_lightcurve', alert_name=alert_names[ii+1]), 
                                alert_names=alert_names,
                                image=fig)
    # Last page.
    elif ii == n_lc - 1:
        return render_template('show_lightcurve_last.html', 
                                home=url_for('query_db'),
                                prev_page=url_for('plot_lightcurve', alert_name=alert_names[ii-1]), 
                                alert_names=alert_names,
                                image=fig)
    # Middle pages.
    else:
        return render_template('show_lightcurve.html', 
                                home=url_for('query_db'),
                                next_page=url_for('plot_lightcurve', alert_name=alert_names[ii+1]), 
                                prev_page=url_for('plot_lightcurve', alert_name=alert_names[ii-1]), 
                                qmax=n_lc-1,
                                alert_names=alert_names,
                                image=fig)

    
def create_figure(time, mag, mag_err, alert_name):
    """
    Plot a lightcurve.
    
    Parameters
    ----------
    time : array-like
        Time (HJD - 2450000)
        
    mag : array-like
        I-band magnitude
    
    mag_err : array-like
        Magnitude uncertainties
        
    alert_name : string
        Name of the alert
        
    Return
    ------
    pngImageB64String : FIXME what is this really?
        This format was chosen so you can pass it into <img src = ...
    """
    #####
    # Figure out limits for plotting the y-axis (magnitude).
    ####
    # MOA alert data is very noisy. We will take the minimum and maximum
    # magnitude range of the observations. But we only use the observations
    # that have error bars in 95% or lower. (Could tweak, I arbitrarily  chose
    # this number to cut out as much junky stuff as possible, but hopefully not 
    # actual data or the peak of the lightcurve.)
    big_err = np.quantile(mag_err, 0.95)
    idx = np.where(mag_err < big_err)[0]
    
    # Get our min and max magnitudes from the less noisy data if necessary.
    if len(idx) == 0:
        # This means that mag_err = big_err. So keep everything.
        ymin = np.min(mag)
        ymax = np.max(mag)
    else:
        ymin = np.min(mag[idx])
        ymax = np.max(mag[idx])
    
    # Change opacity of points depending how many there are.
    npoints = len(time)
    if (npoints <= 5000):
        alpha=0.6
    elif (npoints <= 10000) & (npoints > 5000):
        alpha=0.4
    elif (npoints <= 30000) & (npoints > 10000):
        alpha=0.2
    elif (npoints <= 50000) & (npoints > 30000):
        alpha = 0.05
    else:
        alpha=0.01
    
    # Set up the figure and plot the lightcurve.
    fig = Figure(figsize=(12,6))
    axis = fig.add_subplot(1, 1, 1)
    axis.set_ylim(ymin - 0.2, ymax + 0.2)
    axis.invert_yaxis()
    axis.set_xlabel('HJD - 2450000')
    axis.set_ylabel('I mag')
    axis.errorbar(time, mag, yerr=mag_err, ls='none', marker='.', alpha=alpha, color='k')
    axis.set_title(alert_name)
    
    # Fancy saving stuff: https://stackoverflow.com/questions/61398636/python-flask-matplotlib
    pngImage = io.BytesIO()
    FigureCanvas(fig).print_png(pngImage)
    
    # Encode PNG image to base64 string
    pngImageB64String = "data:image/png;base64,"
    pngImageB64String += base64.b64encode(pngImage.getvalue()).decode('utf8')

    return pngImageB64String


if __name__ == '__main__':
    app.run(port=8000, debug = True)
