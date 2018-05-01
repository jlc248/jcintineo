import logging
from logging.config import dictConfig

def logging_def(logfile=None,MBlimit=10):
  #Creates a stream of output to console by default.
  #With logfile set, will also send the output to that file.

  logging_config = dict(
      version = 1,
      formatters = {
          'f': {'format':'%(asctime)s | %(levelname)s | %(message)s',
                'datefmt':'%Y-%m-%d %H:%M:%S UTC'
          }
      },
      handlers = {
         'h': {'class': 'logging.StreamHandler',
                'formatter': 'f',
                'level': logging.DEBUG
          }
      },
      root = {
          'handlers': ['h'],
          'level': logging.DEBUG,
      },
  )

  if(logfile):
    logging_config['handlers']['file'] = {'class': 'logging.handlers.RotatingFileHandler',
                                          'formatter': 'f',
                                          'filename':logfile,
                                          'maxBytes':MBlimit*1024*1024,
                                          'mode':'a', #use 'a' for append
                                          'backupCount':3
                                         }
    (logging_config['root']['handlers']).append('file')


  return logging_config

