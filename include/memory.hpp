// UjoImro, 2013
// Experimental Research Code for the CARP Project

#ifndef __CARP__MEMORY__HPP__
#define __CARP__MEMORY__HPP__

#include <stdexcept>

namespace carp {
    namespace memory {
        class allocator {
        private:
            int64_t m_pos;
            int64_t m_size;
            
        public:
            static const int64_t memsize = 128;        
            static const int64_t granularity = memsize / 8;

            allocator( int64_t size )
                : m_pos(0)
                , m_size( size )
            {}
                        
            int64_t allocate( int64_t size ) {
                int64_t diff = 0;
                if (m_pos % granularity != 0)
                    diff = granularity - m_pos % granularity;
                m_pos += diff;
                int64_t result = m_pos;
                if (size + m_pos > m_size) throw std::runtime_error("INSUFFICIENT_MEMORY");
                m_pos += size;
                return result;
            }
        };
    }
}

#endif
