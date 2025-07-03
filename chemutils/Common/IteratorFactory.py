import sys
import os
import tempfile
from openeye.oechem import oemolistream, oemolostream, OEWriteMolecule
from chemutils.Common.Util import stdOpen, isStdFile
from chemutils.Common.Const import STD_MOL_EXT


class IteratorFactory:
    """Abstract base class for all iterator factories.
    Given a collection of objects in some form, instances of
    this class should be able to produce (as a factory)
    fresh iterators over the collection.

    Of particular importance is that
    multiple iterator instances should be able to exist
    simultaneously and at different positions to enable
    nested looping over the same collection.

    It effectively encapsulates a streaming data source
    (such as a file) as if it were an in-memory list
    that you could keep retrieving fresh iterators off of.
    """

    def __iter__(self):
        """Primary abstract method where, that returns
        an iterator useable in a "for item in iterator:" construct.
        Based on the "iterable" interface, so no explicit function
        call is needed.  If you want to though, you could do something like

        factory = FileFactory(filename);
        for item in iter(factory):
            print item;

        or just

        for item in factory:
            print item;
        """
        raise NotImplementedError()


class FileFactory(IteratorFactory):
    """Concrete implementation of an IteratorFactory that
    can produce a fresh file cursor iterator over a file
    multiple times.  This requires a reference to the
    filename so that a fresh "open" operation can be called
    on it each time.

    If the filename represents stdin ("-"), then this approach
    will not work and requires a temporary file to be created.
    """

    """File descriptor of temp file created"""
    fd = None
    """Name of the file to open iterators for"""
    filename = None
    """Flag indicating whether a temp file was used"""
    usedTempFile = False

    def __init__(self, aFile):
        """Constructor, taking the filename or file object
        to create iterators for.  If a string filename is given,
        then just use that directly.  Otherwise, if a
        file object is given, or the filename specifies stdin,
        then a temporary file copy will be created.
        """

        if isinstance(aFile, str):
            if not isStdFile(aFile):
                self.filename = aFile
            else:
                aFile = stdOpen(aFile)  # Read stdin as a file object

        if self.filename == None:
            # A filename must not have been passed in (or it was stdin),
            #   create a temporary file to read from then.
            (self.fd, self.filename) = tempfile.mkstemp()
            tempFile = open(self.filename, "w")
            for line in aFile:
                tempFile.write(line)
            tempFile.close()

            self.usedTempFile = True

    def __del__(self):
        """Destructor.  If a temp file was created, then delete it."""
        if self.usedTempFile:
            os.close(self.fd)
            os.remove(self.filename)

    def __iter__(self):
        return open(self.filename)


class oemolistreamFactory(IteratorFactory):
    """Concrete implementation of an IteratorFactory that
    can produce a fresh cursor iterator over an oemolistream.
    This requires a reference to the filename so that a fresh
    "open" operation can be called on it each time.

    If the filename represents stdin ("-"), then this approach
    will not work and requires a temporary file to be created.
    """

    """Pattern for temp filename created"""
    TEMP_NAME = "%d%s"

    """File descriptor and name of the file to open iterators for"""
    fd = None
    filename = None
    """Integer specifying what format the temp file is saved in.
    For example, OEFormat_SDF.  Negative value is a sentinel value
    indicating no temp file is used."""
    tempFileFormat = -1

    def __init__(self, aFile):
        """Constructor, taking the filename or an oemolistream object
        to create iterators for.  If a string filename is given,
        then just use that directly.  Otherwise, if an
        oemolistream object is given, or the filename specifies stdin,
        then a temporary file copy will be created.
        """

        if isinstance(aFile, str):
            if not isStdFile(aFile):
                self.filename = aFile
            else:
                aFile = oemolistream(aFile)  # Read stdin as a file object

        if self.filename == None:
            # A filename must not have been passed in (or it was stdin),
            #   create a temporary file to read from then.
            self.tempFileFormat = aFile.GetFormat()

            # Find the file extension that matches this format code
            fileExt = STD_MOL_EXT[self.tempFileFormat]

            (self.fd, self.filename) = tempfile.mkstemp(suffix=fileExt)
            tempFile = oemolostream(self.filename)
            # tempFile.SetFormat(self.tempFileFormat)

            for mol in aFile.GetOEGraphMols():
                OEWriteMolecule(tempFile, mol)
            tempFile.close()

    def __del__(self):
        """Destructor.  If a temp file was created, then delete it."""
        if self.tempFileFormat > -1:
            os.close(self.fd)
            os.remove(self.filename)

    def __iter__(self):
        oeis = oemolistream(self.filename)
        # if self.tempFileFormat > -1:
        #    oeis.SetFormat(self.tempFileFormat)  # Need to specify a molecule format if temp file used
        return oeis.GetOEGraphMols()


def main(argv):
    """Main method, callable from command line"""
    # print("This is an abstract base class and should not be invoked directly.")
    # print("Use one of the concrete Kernel classes instead.")
    sys.exit(-1)


if __name__ == "__main__":
    main(sys.argv)
