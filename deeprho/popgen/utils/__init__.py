import os
import jpype
import jpype.imports
from jpype.types import *
from ... import popgen
from .treeutils import *
from .utils import *
from .readers import *
from .simulator import *
from .statistics import *
package_dirname = os.path.dirname(popgen.__file__)

# start JVM to import .jar libs
if not jpype.isJVMStarted():
    jpype.startJVM(classpath=[os.path.join(package_dirname, 'libs/*')])

from .javautils import *




