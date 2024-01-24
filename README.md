Microlensing Alert Database (MAD)

Initially created by Casey Lam, further developed by Anette Brecko and Macy Huston

# MAD 
Microlensing Alerts Database (MAD) is a Flask webapp that allows the user to 
interact with a database of event parameters and lightcurves from the
[MOA](http://www.massey.ac.nz/~iabond/moa/alerts/), 
[OGLE](https://ogle.astrouw.edu.pl/ogle4/ews/), and
[KMTNet](https://kmtnet.kasi.re.kr/~ulens/)
microlensing alert websites.

With MAD, the user can query the database and download the results as a CSV file.

In addition, the user can query the database to choose lightcurves to view.

### Background and motivation:
Our research group is looking for isolated stellar mass black holes with microlensing.
We want to choose promising black hole candidates from photometric microlensing 
alert systems to follow up astrometrically for our HST Cycle 29 program.
The way we've done this in the past is to look through individual alert pages for events of interest.
However, this is not ideal for many reasons:
1. It is very time consuming to click through each alert page one at a time to look at lightcurves.
Although all the fit parameters are located on the home page, MOA and OGLE report uncertainties 
on their parameters which are only shown on the individual event pages.
2. The active alert page for the current year are updated on a regular (nightly) basis.
This means that the reported best-fit microlensing model parameters also get updated, 
so an event that looks promising might not be so good two weeks later or vice versa.  
3. Sometimes we want to re-analyze a set of lightcurves, 
and having to go page by page to download what we want is also time consuming.
4. Combining all the issues above leads to a general mess.
This also makes the process likely less thorough.
5. For MOA in particular, the lightcurves provided for visualization on the alert website
are shown as difference image lightcurves, i.e. delta flux vs. time. 
However, this is impossible to interpret in the sense that we care about (i.e. what is the
difference between the baseline and peak magnification?)

MAD was created to help alleviate these issues.
By aggregating all the alerts into a database, we can easily query a table of 
alert parameters to find the events that are most promising.
Then, we can use JOIN queries to grab the photometry data for those events
to re-analyze and visualize ourselves.
The database can then also be easily updated to grab the latest changes, and
the same query can be re-run, to see what changes.

### File structure of MAD
There are three python files:
1. `app.py`. 
This is the Flask app itself.
2. `query_alerts.py`. 
This is a module that contains a bunch of functions 
that can be used to  populate the database.
It primarily consists of functions that scrape microlensing event alerts 
and lightcurves from the OGLE, KMTNet, and MOA websites.
3. `populate_database.py`.
This is a very short script that calls `query_alerts.py` that populates the database.
This can be edited by the user depending on what set of alerts and photometry
they are particularly interested in.
4. `fitting_utils.py`.
This module is called from app.py when the output database
is downloaded as a JSON file, and it downloads the lightcurve of each event
5. `mad_munge.py`
Sets priors for model fitting
7. templates
Folder with all of the files for HTML display pages
8. `run.py` in bagle_fits folder
Takes the JSON file from the output database and feeds it as in input
into Bayesian Analysis of Gravitational Lensing Events (BAGLE)

### Workflow
Note that before the Flask app (`app.py`) can be run, `populate_database.py` must have been run.

By running `populate_database.py`, a database called `microlensing.db` is created.
It consists of up to two tables, `alerts` and `photometry`.
The primary key for the `alerts` table is `alert_name`. 
The `alerts` table also contains the microlensing parameters reported by the alert system,
as well as URLs to the original alert pages, among other information.
In the `photometry` table, `alert_name` is the foreign key.
Each row in the `photometry` table also contains the date, magnitude, magnitude uncertainty, 
and telescope information for each observation.

### Future work
1. Figure out why populate_database currently only runs in ipython, not command line.
2. Parallelize lightcurve download
3. Automating database updates via CronJob
