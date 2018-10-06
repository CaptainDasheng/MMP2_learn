'''
Created on Nov 5, 2017

@author: Yuhao Fu
'''

class LogHandler():
    """
    handler of log information.
    """
    def __init__(self, path=None, **kwargs):
        """
        Arguments:
            path: path of log. i.e. /home/xx/xx/run.log
            
            kwargs:
                isCovered (default=False): whether to cover the old log information.
        """
        import os
        import warnings
        
        if path is None:
            self.path='run.log'
        else:
            self.path=path
            if os.path.basename(path).lower().startswith('.'):
                warnings.warn("log file is a hidden file: '%s'" %os.path.basename(path))
            if not os.path.basename(path).lower().endswith('.log'):
                warnings.warn('suffix of log file is not .log')
        
        # whether to clear the old log information.
        isCovered=False
        if 'isCovered' in kwargs:
            isCovered=kwargs['isCovered']
            if not isinstance(isCovered, bool):
                raise ValueError('isCovered is not belong to boolean value')
        if isCovered:
            os.remove(self.path)
            
    def write(self, string):
        """
        write a message to log file.
        
        Arguments:
            string: log message.
        """
        from utils.fetch import get_time
        
        output=open(self.path, 'a')
        output.write('%s %s\n' %(get_time(), string))
        output.close()