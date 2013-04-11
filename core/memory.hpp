// UjoImro, 2013
// Experimental Research Code for the CARP Project

#ifndef __CARP__MEMORY__HPP__
#define __CARP__MEMORY__HPP__

#include <bitset>
#include <memory>
#include <cstddef>
#include <utility>

namespace carp {

    template <class T0>
    class memory {
    private:

        typedef std::bitset<8*sizeof(int64_t)> position_t;
                
        int64_t m_size;        

        // three statuses:
        // free - the block is not allocated
        // broken - the block is broken into subblocks (left and right are not NULL)
        // allocated - the block is used for allocation
        enum status_t { free, broken, allocated };

        int64_t
        log2ceil( const int64_t & val ) {
            int64_t result;
            for ( int64_t = 0; (1ul<<result) < val; result++);
            return result;
        } // log2ceil 
        
        class block_t {
        public:
            status_t m_status;
            int64_t  m_level;
            std::shared_ptr<block_t> m_left;
            std::shared_ptr<block_t> m_right;
            
            block_t( const int64_t & level ) : m_status(free), m_level(level) { assert(level>=0); };

            std::pair<bool, int64_t>
            acquire( const int64_t & level ) {
                if ( m_level < level )
                    return make_pair(false, 0);
                
                if ( m_status == allocated ) 
                    return make_pair(false, 0);

                if ( m_status == free ) {
                    if ( m_level == level ) { // we are allocating this block for the request
                        m_status = allocated;
                        return std::make_pair( true, 0 );
                    }
                    else { // NOT m_level == level (m_level > level ) 
                        assert(m_level>level);                    
                        split();
                    }
                } // m_status == free

                // if we have got here then m_status == broken AND m_level>level (so m_level>0)
                assert(m_level > level);
                auto left_node = left->acquire(level);
                if (left_node.first)
                    return std::make_pair( true, left_node.second );

                auto right_node = right->acquire(level);
                if (right_node.first)
                    return std::make_pair( true, (1<<level) + right_node.second );

                return std::make_paid(false, 0);
            } // acquire

            bool split() {
                assert( m_level > 0 );
                assert( m_status == free );
                m_status = broken;
                m_left.reset( new block(m_level-1) );
                m_right.reset( new block(m_level-1) );
            } // split
                        
            void merge() {
                if ( ( m_status == broken ) and
                     ( left.status  == free ) and
                     ( right.status == free ) ) {
                    left.reset();
                    right.reset();
                    m_status = free;
                } // if
            } // merge

            void release( const position_t & pos ) throw ( memory::exception& ){
                if (m_status == free) throw memory::exception(FREEING_UNALLOCATED_MEMORY);
                    
                if (m_status == allocated) {
                    m_status == free;
                    return;                    
                }

                // here the status is broken
                assert(m_status==broken);
                
                if (pos.test(m_level))
                    left->release();
                else // NOT nextbig
                    right->release();

                merge(); // we merge the blocks if it's possible                

                return;                
            } // release
                        
        }; // class block

        block pool;        
        
    public:
        
        typedef T0 value_type;        
        static const int64_t memsize = 128;        
        static const int64_t granularity = memsize / 8;
        
        enum error_t { INSUFFICIENT_MEMORY, FREEING_UNALLOCATED_MEMORY };
        
        class exception {
            public error_t error;
        }; // class exception 

        memory ( const int64_t & numel )
            : m_size( sizeof(value_type) * size ),
              pool( log2ceil(m_size) ) { } // memory

        int64_t
        allocate( const int64_t & numel ) throw( memory::exception& ) {
            int64_t size = numel * sizeof(value_type);            

            level = log2ceil( size / memsize );
            auto node = pool.allocate(level);
            if (!node.first) throw memory::exception(INSUFFICIENT_MEMORY);

            return 128 * node.second / sizeof(value_type);
        } // allocate

        void
        release( const int64_t & elpos ) throw ( memory::exception& ) {
            // int64_t size = elpos * sizeof(value_type);

            assert( elpos * sizeof(value_type) % 128 == 0 );
            
            position_t position( elpos * sizeof(value_type) / 128 );
            pool.release(position);            

            return;
        } // release
        
    }; // memory
        
} // namespace CARP






#endif /* __CARP__MEMORY__HPP__ */
// LuM end of file
