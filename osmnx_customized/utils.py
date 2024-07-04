"""General utility functions."""

import datetime as dt
import logging as lg
import os
import sys
import unicodedata
import warnings
from contextlib import redirect_stdout
from pathlib import Path

from . import settings


def citation():
    """
    Print the OSMnx package's citation information.

    Boeing, G. 2017. OSMnx: New Methods for Acquiring, Constructing, Analyzing,
    and Visualizing Complex Street Networks. Computers, Environment and Urban
    Systems, 65, 126-139. https://doi.org/10.1016/j.compenvurbsys.2017.05.004

    Returns
    -------
    None
    """
    cite = (
        "Citation:\n\n"
        "Boeing, G. 2017. OSMnx: New Methods for Acquiring, "
        "Constructing, Analyzing, and Visualizing Complex Street "
        "Networks. Computers, Environment and Urban Systems, 65, 126-139. "
        "https://doi.org/10.1016/j.compenvurbsys.2017.05.004\n\n"
        "BibTeX entry for LaTeX users:\n\n"
        "@article{boeing_osmnx_2017,\n"
        "    title = {{OSMnx}: {New} {Methods} for {Acquiring}, "
        "{Constructing}, {Analyzing}, and {Visualizing} {Complex} "
        "{Street} {Networks}},\n"
        "    volume = {65},\n"
        "    doi = {10.1016/j.compenvurbsys.2017.05.004},\n"
        "    number = {126-139},\n"
        "    journal = {Computers, Environment and Urban Systems},\n"
        "    author = {Boeing, Geoff},\n"
        "    year = {2017}\n"
        "}"
    )

    print(cite)


def ts(style="datetime", template=None):
    """
    Get current timestamp as string.

    Parameters
    ----------
    style : string {"datetime", "date", "time"}
        format the timestamp with this built-in template
    template : string
        if not None, format the timestamp with this template instead of one of
        the built-in styles

    Returns
    -------
    ts : string
        the string timestamp
    """
    if template is None:
        if style == "datetime":
            template = "{:%Y-%m-%d %H:%M:%S}"
        elif style == "date":
            template = "{:%Y-%m-%d}"
        elif style == "time":
            template = "{:%H:%M:%S}"
        else:  # pragma: no cover
            raise ValueError(f'unrecognized timestamp style "{style}"')

    ts = template.format(dt.datetime.now())
    return ts


def config(
    all_oneway=settings.all_oneway,
    bidirectional_network_types=settings.bidirectional_network_types,
    cache_folder=settings.cache_folder,
    cache_only_mode=settings.cache_only_mode,
    data_folder=settings.data_folder,
    default_accept_language=settings.default_accept_language,
    default_access=settings.default_access,
    default_crs=settings.default_crs,
    default_referer=settings.default_referer,
    default_user_agent=settings.default_user_agent,
    imgs_folder=settings.imgs_folder,
    log_console=settings.log_console,
    log_file=settings.log_file,
    log_filename=settings.log_filename,
    log_level=settings.log_level,
    log_name=settings.log_name,
    logs_folder=settings.logs_folder,
    max_query_area_size=settings.max_query_area_size,
    memory=settings.memory,
    nominatim_endpoint=settings.nominatim_endpoint,
    nominatim_key=settings.nominatim_key,
    osm_xml_node_attrs=settings.osm_xml_node_attrs,
    osm_xml_node_tags=settings.osm_xml_node_tags,
    osm_xml_way_attrs=settings.osm_xml_way_attrs,
    osm_xml_way_tags=settings.osm_xml_way_tags,
    overpass_endpoint=settings.overpass_endpoint,
    overpass_rate_limit=settings.overpass_rate_limit,
    overpass_settings=settings.overpass_settings,
    requests_kwargs=settings.requests_kwargs,
    timeout=settings.timeout,
    use_cache=settings.use_cache,
    useful_tags_node=settings.useful_tags_node,
    useful_tags_way=settings.useful_tags_way,
):
    """
    Do not use: deprecated. Use the settings module directly.

    Parameters
    ----------
    all_oneway : bool
        deprecated
    bidirectional_network_types : list
        deprecated
    cache_folder : string or pathlib.Path
        deprecated
    data_folder : string or pathlib.Path
        deprecated
    cache_only_mode : bool
        deprecated
    default_accept_language : string
        deprecated
    default_access : string
        deprecated
    default_crs : string
        deprecated
    default_referer : string
        deprecated
    default_user_agent : string
        deprecated
    imgs_folder : string or pathlib.Path
        deprecated
    log_file : bool
        deprecated
    log_filename : string
        deprecated
    log_console : bool
        deprecated
    log_level : int
        deprecated
    log_name : string
        deprecated
    logs_folder : string or pathlib.Path
        deprecated
    max_query_area_size : int
        deprecated
    memory : int
        deprecated
    nominatim_endpoint : string
        deprecated
    nominatim_key : string
        deprecated
    osm_xml_node_attrs : list
        deprecated
    osm_xml_node_tags : list
        deprecated
    osm_xml_way_attrs : list
        deprecated
    osm_xml_way_tags : list
        deprecated
    overpass_endpoint : string
        deprecated
    overpass_rate_limit : bool
        deprecated
    overpass_settings : string
        deprecated
    requests_kwargs : dict
        deprecated
    timeout : int
        deprecated
    use_cache : bool
        deprecated
    useful_tags_node : list
        deprecated
    useful_tags_way : list
        deprecated

    Returns
    -------
    None
    """
    warnings.warn(
        "The `utils.config` function is deprecated and will be removed in a "
        "future release. Instead, use the `settings` module directly to "
        "configure a global setting's value. For example, "
        "`ox.settings.log_console=True`."
    )

    # set each global setting to the argument value
    settings.all_oneway = all_oneway
    settings.bidirectional_network_types = bidirectional_network_types
    settings.cache_folder = cache_folder
    settings.cache_only_mode = cache_only_mode
    settings.data_folder = data_folder
    settings.default_accept_language = default_accept_language
    settings.default_access = default_access
    settings.default_crs = default_crs
    settings.default_referer = default_referer
    settings.default_user_agent = default_user_agent
    settings.imgs_folder = imgs_folder
    settings.log_console = log_console
    settings.log_file = log_file
    settings.log_filename = log_filename
    settings.log_level = log_level
    settings.log_name = log_name
    settings.logs_folder = logs_folder
    settings.max_query_area_size = max_query_area_size
    settings.memory = memory
    settings.nominatim_endpoint = nominatim_endpoint
    settings.nominatim_key = nominatim_key
    settings.osm_xml_node_attrs = osm_xml_node_attrs
    settings.osm_xml_node_tags = osm_xml_node_tags
    settings.osm_xml_way_attrs = osm_xml_way_attrs
    settings.osm_xml_way_tags = osm_xml_way_tags
    settings.overpass_endpoint = overpass_endpoint
    settings.overpass_rate_limit = overpass_rate_limit
    settings.overpass_settings = overpass_settings
    settings.timeout = timeout
    settings.use_cache = use_cache
    settings.useful_tags_node = useful_tags_node
    settings.useful_tags_way = useful_tags_way
    settings.requests_kwargs = requests_kwargs


def log(message, level=None, name=None, filename=None):
    """
    Write a message to the logger.

    This logs to file and/or prints to the console (terminal), depending on
    the current configuration of settings.log_file and settings.log_console.

    Parameters
    ----------
    message : string
        the message to log
    level : int
        one of Python's logger.level constants
    name : string
        name of the logger
    filename : string
        name of the log file, without file extension

    Returns
    -------
    None
    """
    if level is None:
        level = settings.log_level
    if name is None:
        name = settings.log_name
    if filename is None:
        filename = settings.log_filename

    # if logging to file is turned on
    if settings.log_file:
        # get the current logger (or create a new one, if none), then log
        # message at requested level
        logger = _get_logger(level=level, name=name, filename=filename)
        if level == lg.DEBUG:
            logger.debug(message)
        elif level == lg.INFO:
            logger.info(message)
        elif level == lg.WARNING:
            logger.warning(message)
        elif level == lg.ERROR:
            logger.error(message)

    # if logging to console (terminal window) is turned on
    if settings.log_console:
        # prepend timestamp
        message = f"{ts()} {message}"

        # convert to ascii so it doesn't break windows terminals
        message = (
            unicodedata.normalize("NFKD", str(message)).encode("ascii", errors="replace").decode()
        )

        # print explicitly to terminal in case jupyter notebook is the stdout
        if getattr(sys.stdout, "_original_stdstream_copy", None) is not None:
            # redirect captured pipe back to original
            os.dup2(sys.stdout._original_stdstream_copy, sys.__stdout__.fileno())
            sys.stdout._original_stdstream_copy = None
        with redirect_stdout(sys.__stdout__):
            print(message, file=sys.__stdout__, flush=True)


def _get_logger(level, name, filename):
    """
    Create a logger or return the current one if already instantiated.

    Parameters
    ----------
    level : int
        one of Python's logger.level constants
    name : string
        name of the logger
    filename : string
        name of the log file, without file extension

    Returns
    -------
    logger : logging.logger
    """
    logger = lg.getLogger(name)

    # if a logger with this name is not already set up
    if not getattr(logger, "handler_set", None):

        # get today's date and construct a log filename
        log_filename = Path(settings.logs_folder) / f'{filename}_{ts(style="date")}.log'

        # if the logs folder does not already exist, create it
        log_filename.parent.mkdir(parents=True, exist_ok=True)

        # create file handler and log formatter and set them up
        handler = lg.FileHandler(log_filename, encoding="utf-8")
        formatter = lg.Formatter("%(asctime)s %(levelname)s %(name)s %(message)s")
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.setLevel(level)
        logger.handler_set = True

    return logger
