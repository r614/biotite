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
        self._content = []
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
    
    @classmethod
    def write(self, file):
        """
        Write contents into a BinaryCIF file.
        
        Parameters
        ----------
        file : file-like object or str
            The file to be written to.
            Alternatively, a file path can be supplied.
        """

        serializedObject = self._serialize()
        
        packed_bytes = msgpack.packb(
            serializedObject, use_bin_type=True, default=_encode_numpy
        )
        if isinstance(file, str):
            with open(file, "wb") as f:
                f.write(packed_bytes)
        else:
            if not is_binary(file):
                raise TypeError("A file opened in 'binary' mode is required")
            file.write(packed_bytes)

    def _serialize(self):
        res = []
        for block in self._content: 
            res.append(block._serialize())
        return res

    def _parse(self, raw_blocks):
        dataBlocks = []

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
        Category: A Category object. 
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
            self._block_map(block).set_category(category_dict)
        else:
            block = "data_" + category_dict.get("name")
            new_block = DataBlock(block)
            new_block.set_category(category_dict)
            self._content.append(new_block)
            self._block_map[block] = new_block

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
    
    def get_category(self, category_name):
        return self._category_map.get(category_name)

    def set_category(self, category_name, category_dict):
        if self.get_category(category_name):
            self.get_category(category_name).set_category(category_dict)
        
        new_category = Category(category_dict['name'],
                                category_dict['rowCount'],
                                category_dict['columns'])
        self._categories.append(new_category)
        self._category_map[category_dict['name']] = new_category

    def _serialize(self):
        res = []
        for category in self._categories:
            res.append(category._serialize())


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

    def _serialize(self):
        res = []
        for column in columns:
            res.append(columns._serialize())
        return res

    def set_category(self, category_dict):
        self._columns = self._build_columns(category_dict.get('columns'))
        self._name = category_dict.get('name')
        self._row_count = category_dict.get('rowCount')

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
    
    def _serialize(self):
        if not mask:
            return {
                'data': self._data._serialize()
            }
        return {
            'data': self._data._serialize(),
            'mask': self._mask._serialize()
        }


class Data(MutableMapping):
    def __init__(self, data=None, encoding=None, decoded=False):
        self._data = data
        self._encoding = encoding
        # self.decoded = self._decode(self._data)
        self.is_decoded = decoded
        self.decoded_data = None 

    def _serialize(self):
        if not self.is_decoded:
            return {
                'data': self.data
                'encoding': self.encoding
                }
        return self.encode(self.data, self._encoding)

    def __decode(self):
        // TODO

        self.decoded_data = decode_data({'data': self._data,
                                        'encoding': self._encoding})
        return None
    
    def _encode(self):
        // TODO 
        return None



