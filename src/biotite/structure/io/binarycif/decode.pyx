# This source code is part of the Biotite package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

__name__ = "biotite.structure.io.binarycif"
__author__ = "Roshan Pawar"
__all__ = ["decode", "Decoder"]

cimport cython
cimport numpy as np

import numpy as np
# Maybe writing these functions for BinaryCIF again won't be a bad idea because of a future change in spec?
#from .decode import _decode_run_length as common_decode_run_length 
#from .decode import _decode_packed as common_decode_packed 


ctypedef np.int8_t int8
ctypedef np.int16_t int16
ctypedef np.int32_t int32
ctypedef np.uint8_t uint8
ctypedef np.uint16_t uint16
ctypedef np.uint32_t uint32
ctypedef np.uint64_t uint64
ctypedef np.float32_t float32
ctypedef np.float64_t float64


ctypedef fused deltaType:
        int8
        int16
        int32


def decode_data(dict Data):
    cdef uint8[:] data_array = Data['data']
    for encoder in Data['encoding'][::-1]:
        if 'kind' not in encoder:
            raise ValueError('Error parsing encoder, "kind" field not found.')
        binary_data = Decoder.decodeStep(data_array, encoder)
    return data_array

class Decoder:

    def decodeStep(self, uint8[:] data, dict encoder):
            
        if encoder['kind'] == 'StringArray': return self._stringArray(data, encoder)
        elif encoder['kind'] == 'IntervalQuantization': return self._intervalQuantization(data, encoder)
        elif encoder['kind'] == 'ByteArray': return self._resolveByteArray(data, encoder) 
        elif encoder['kind'] == 'Delta': return self._delta(data, encoder)
        #elif encoder['kind'] == 'RunLength': return common_decode_run_length(data)
        #elif encoder['kind'] == 'FixedPoint': return common_decode_packed(data)
        elif encoder['kind'] == 'IntegerPacking': return self._integerPacking(data, encoder)
        else: return data

    def _stringArray(self, uint8[:] data, dict encoding):
        cdef str nested_data = encoding['stringData']
        cdef np.ndarray offsets = decode_data({ 'encoding': encoding['offsetEncoding'], 'data': encoding['offsets']  })
        cdef np.ndarray indices = decode_data({ 'encoding': encoding['dataEncoding'], 'data': data  })
        cache = list()
        result = [None for _ in indices]
        cdef int offset = 0
        for i in indices:
            if (i < 0):
                offset += 1
                result[offset] = None
                continue
            v = cache[i]
            if v is None:
                v = nested_data[offsets[i]:offsets[i+1]]
                cache[i] = v
        return result 
    
    def _getFloatArray(self, dict encoding, int size):
        if encoding['srcType'] == 'Float32':
            return np.ndarray(size, dtype=np.float32)
        return np.ndarray(size, dtype=np.float64)

    def _resolveIntType(self, uint8[:] data, str src):
        if src == 'Int8': return np.int8
        elif src == 'Int16': return np.int16
        elif src == 'Int32': return np.int32
        elif src == 'Int64': return np.int64
        elif src == 'Uint8': return np.uint8
        elif src == 'Uint16': return np.uint16
        elif src == 'Uint32': return np.uint32
        elif src == 'Uint64': return np.uint64
        else: return None

    def _intervalQuantization(self, int32[:] data, dict encoding):
        cdef int n = len(data)
        output = self._getFloatArray(encoding, n)
        delta = (encoding['max'] - encoding['min'])/(encoding['numSteps'] - 1)
        for i in range(n):
            output[i] = encoding['min'] + delta * data[i] 
        return output

    def _resolveByteArray(self, uint8[:] data, dict encoding):
        inttype = self._resolveIntType(encoding['type'])

        return np.array(data).astype(inttype)
    
    
    def _delta(self, deltaType[:] array, dict encoding):
        int_type = self._resolveIntType(encoding['srcType'])
        return np.cumsum(array, dtype=int_type)
    
    def _fixedPoint(self, int32[:] data, dict encoding):
        cdef int n = len(data)
        output = self._getFloatArray(encoding, n)
        cdef double f = 1/encoding['factor']
        for i in range(n):
            output[i] = f * data[i]
        return output

