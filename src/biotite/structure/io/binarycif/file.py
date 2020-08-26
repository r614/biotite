# This source code is part of the Biotite package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

__name__ = "biotite.structure.io.binarycif"
__author__ = "Roshan Pawar"
__all__ = ["BinaryCIFFile"]

from collections.abc import MutableMapping
import numpy as np
import msgpack
from ....file import TextFile
# from ..encode import *
from .decode import decode_data


class BinaryCIFFile(TextFile, MutableMapping):
    """
    This class represents a BinaryCIF file.
    
    """
    
    def __init__(self):
        super().__init__()
        self._content = {}
        self._content["encoder"] = "UNKNOWN"
        self._block_map = {}
     
    @classmethod
    def read(self, file):
        """
        Read a BCIF file.
        
        Parameters
        ----------
        file : file-like object or str
            The file to be read.
            Alternatively a file path can be supplied.
        
        Returns
        -------
        file_object : BinaryCIFFile
            The parsed file.
        """
        bcif_file = BinaryCIFFile()
        raw_content = bcif_file._content
        # File name
        if isinstance(file, str):
            with open(file, "rb") as f:
                raw_content = msgpack.unpackb(
                    f.read(), use_list=True, raw=False
                )
        # File object
        else:
            if not is_binary(file):
                raise TypeError("A file opened in 'binary' mode is required")
            raw_content = msgpack.unpackb(
                file.read(), use_list=True, raw=False
            )
        
        bcif_file._content = bcif_file._parse(raw_content['dataBlocks'])
        
        return bcif_file

    def _parse(self, raw_blocks):
        dataBlocks = [];

        for block in raw_blocks:
            current = DataBlock(block["header"], block["categories"])
            dataBlocks.append(DataBlock(block["header"], block["categories"]))
            self._block_map[block['header']] = current

        return dataBlocks

    def get_block_headers(self):
        """
        Get the headers of all data blocks in the file.
        
        Returns
        -------
        blocks : list
            List of data block header.
        """
        blocks = set()
        for block in self._content():            
            blocks.add(block._header)
        return sorted(blocks)
    
    def get_category(self, category, block=None):
        """
        Get the dictionary for a given category.
        
        Parameters
        ----------
        category : string
            The name of the category. The leading underscore is omitted.
        block : string, optional
            The header of the data block. Default is the first
            (and most times only) data block of the file.
            
        Returns
        -------
        category_dict : dict of (str or ndarray, dtype=str) or None
            A entry keyed dictionary. The corresponding values are
            strings or array of strings
            Returns None, if the data block does not contain the given
            category.
        """
        if block is None:
            block = self.get_block_headers()[0]

        block = self._block_map.get(block)
        if not block:
            return None

        category = block._category_map.get(category)
        if category is None:
            return None
        return category.__dict__()
    
    def set_category(self, category, category_dict, block=None):
        """
        Set the content of a category.
        
        If the category is already exisiting, all lines corresponding
        to the category are replaced. Otherwise a new category is
        created and the lines are appended at the end of the data block.
        
        Parameters
        ----------
        category : string
            The name of the category. The leading underscore is omitted.
        category_dict : dict
            The category content. The dictionary must have strings
            (subcategories) as keys and strings or :class:`ndarray`
            objects as values.
        block : string, optional
            The name of the data block. Default is the first
            (and most times only) data block of the file. If the
            block is not contained in the file yet, a new block is
            appended at the end of the file.
        """
        if not block:
            block = self.get_block_headers()[0]
        
        if block: 
            self._categories[(category, block)] = category_dict
        else:
            block = "data_" + category_dict.get("name")
            self._categories[(category, block)] = category_dict

    def __delitem__():
        return 

    def __getitem__():
        return 

    def __iter__():
        return

    def __len__():
        return

    def __setitem__():
        return 


class DataBlock(MutableMapping):
    def __init__(self, header=None, categories=[]):
        self._header = header
        self._categories = self._build_categories(categories)
        self._category_map = {}
    
    def _build_categories(self, categories):
        if not categories:
            return categories
        res = []
        for category in categories:
            current = Category(category['name'], category['rowCount'],
                               category['columns'])
            res.append(current)
            self._category_map[category['name']] = current
        return res


class Category(MutableMapping):
    def __init__(self, name, row_count, raw_columns):
        self._columns = self._build_columns(raw_columns)
        self._name = None
        self._row_count = None

    def _build_columns(self, raw_columns):
        res = []
        for col in raw_columns:
            res.append(Column(col['name'], col['data'], col['mask']))
        return res

    def __delitem__():
        return 

    def __getitem__():
        return 

    def __iter__():
        return

    def __len__():
        return

    def __setitem__():
        return 


class Column(MutableMapping):
    def __init__(self, name=None, data=None, mask=None):
        self._data = Data(data['data'], data['encoding'])
        self._name = name
        self._mask = Data(mask.get('data'), mask.get('encoding'))


class Data(MutableMapping):
    def __init__(self, data=None, encoding=None, decoded=False):
        self._data = data
        self._encoding = encoding
        self.decoded = decoded
        self.is_decoded = decoded

    def __decode(self):
        // TODO
        return None
        



