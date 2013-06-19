// UjoImro, 2013
// Experimental Research Code for the CARP Project

#ifndef __CARP__MEMORY__HPP__
#define __CARP__MEMORY__HPP__

#include <bitset>
#include <memory>
// #include <string>
#include <cstddef>
#include <cstdint>
#include <utility>
// #include <iostream>
#include <assert.h>
#include <algorithm>
#include <boost/preprocessor.hpp>

#include "utility.hpp"

// template <class T0>
// void
// print_pos( T0 & pos, std::string name ) {
//     std::cout << name << " = ";
//     for (int q=0; q<8; q++ )
//         std::cout << pos[q] << ", ";

//     std::cout << "..." << std::endl;
    
//     return;
    
// } // print_pos


// #include "errors.hpp"

namespace carp {

    namespace memory {        
    
        enum error_t { INSUFFICIENT_MEMORY, FREEING_UNALLOCATED_MEMORY };
    
        class exception : public std::exception {
        public:
            error_t error;

            exception( const error_t & error ) : error(error) { }
        
            virtual const char* what() const noexcept /* override */ {
                switch (error) {
                case INSUFFICIENT_MEMORY:
                    return "INSUFFICIENT_MEMORY";
                    break;
                case FREEING_UNALLOCATED_MEMORY:
                    return "FREEING_UNALLOCATED_MEMORY";
                    break;
                
                } // switch
                return "";
            
            }
        
        }; // class exception 


        class allocator {
        public:

            static const int64_t memsize = 128;        
            static const int64_t granularity = memsize / 8;
            static const int64_t log2granularity = 4;
            
            class sizer {
            public:
                int64_t numel;
                int64_t size;
                template <class MT0> sizer( int64_t numel, MT0 ) : numel(numel), size(sizeof(MT0)) { }

                int64_t addr() const {
                    return numel * size;                    
                }
                
            }; // class sizer
                
            virtual int64_t grossallocated() const = 0;
            
            virtual int64_t netallocated() const = 0;
            
            int64_t virtual allocate( sizer allocsize ) throw( carp::memory::exception& ) = 0;
            
            void virtual release( sizer elpos ) throw ( carp::memory::exception& ) {
                assert(false);
            } // release

            virtual ~allocator() { };
        
        }; // class allocator
    

        
        class buddy : public allocator {
        private:

            typedef std::bitset<8*sizeof(int64_t)> position_t;
                
            int64_t m_size;
            int64_t m_net_allocated;
            int64_t m_gross_allocated;
        
            // three statuses:
            // free - the block is not allocated
            // broken - the block is broken into subblocks (left and right are not NULL)
            // allocated - the block is used for allocation
            enum status_t { free, broken, allocated };

            int64_t
            log2ceil( const int64_t & val ) {
                int64_t result;
                for ( result = 0; (1ul<<result) < val; result++);
                return result;
            } // log2ceil 
        
            class block_t {
            public:
                status_t m_status;
                int64_t  m_level;
                std::shared_ptr<block_t> m_left;
                std::shared_ptr<block_t> m_right;
            
                block_t( const int64_t & level ) : m_status(free), m_level(level) {
                    // PRINT("block_t::block_t");
                    assert(level>=0);
                    // PRINT(m_level);
                }; // block_t::block_t

                std::pair<bool, int64_t>
                allocate( const int64_t & level ) {
                    // PRINT("block_t::allocate");
                    // PRINT(m_level);
                    // PRINT(level);
                                   
                    if ( m_level < level )
                        return std::make_pair(false, 0);
                
                    if ( m_status == allocated ) 
                        return std::make_pair(false, 0);

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

                    assert( m_status == broken );

                    if (m_level==level)
                        // the block is of correct size, but it is split
                        // up (broken up), which means, that at least half
                        // of it is allocated, which means that the block
                        // cannot be used to accomodate a memory of the
                        // requested size
                        return std::make_pair( false, 0 );
                
                
                    // if we have got here then m_status == broken AND m_level>level (so m_level>0)
                    assert(m_level > level);
                    auto left_node = m_left->allocate(level);
                    // PRINT(left_node.first);
                    // PRINT(left_node.second);
                 
                    if (left_node.first)
                        return std::make_pair( true, left_node.second );
                
                    auto right_node = m_right->allocate(level);
                    // PRINT(right_node.first);
                    // PRINT(right_node.second);
                 
                    if (right_node.first)
                        return std::make_pair( true, (1<<(m_level-1)) + right_node.second );
                
                    return std::make_pair(false, 0);
                } // allocate

                void split() {
                    // PRINT("block_t::split");
                    assert( m_level > 0 );
                    assert( m_status == free );
                    m_status = broken;
                    m_left.reset( new block_t(m_level-1) );
                    m_right.reset( new block_t(m_level-1) );
                } // split
                        
                void merge() {
                    // PRINT("block_t::merge");
                
                    if ( ( m_status == broken ) &&
                         ( m_left->m_status  == free ) &&
                         ( m_right->m_status == free ) ) {
                        m_left.reset();
                        m_right.reset();
                        m_status = free;
//                    PRINT("merge occurred");
                    
                    } // if
                } // merge

                int64_t
                release( const position_t & pos) throw ( carp::memory::exception& ){
                    if (m_status == free) throw carp::memory::exception(FREEING_UNALLOCATED_MEMORY);
                    
                    if (m_status == allocated) {
//                    PRINT("found an allocated node");
                        m_status = free;
                        return m_level;
                    }

                    // here the status is broken
                    assert(m_status==broken);

//                PRINT(m_level);

                    auto released_size = 0;
                
                    if (! pos.test(m_level-1)) {
//                    PRINT("went left");
                    
                        released_size = m_left->release(pos);
                    }
                    else {// NOT nextbig
//                    PRINT("went right");
                        released_size = m_right->release(pos);
                    }
                
                    merge(); // we merge the blocks if it's possible                

                    return released_size;
                } // release
                        
            }; // class block_t

            block_t pool;        
        
        public:
        
            buddy ( sizer allocsize )
                : m_size( allocsize.addr() ),
                  m_net_allocated(0),
                  m_gross_allocated(0),
                  pool( log2ceil( allocsize.addr() ) - log2granularity ) {
                // only sizes which divide 128bits are supported
                // PRINT("memory::memory");
                assert( (granularity % allocsize.size) == 0);
//            PRINT(numel);
            
//            PRINT( log2ceil( numel * sizeof(value_type) ) - log2granularity );
            } // memory

            int64_t
            grossallocated() /* override */ const {
                return m_gross_allocated;
            }

            int64_t netallocated() /* override */ const {
                return m_net_allocated;            
            }

            int64_t
            allocate( sizer allocpos ) throw( carp::memory::exception& ) /* override */ {
                assert(allocpos.numel>0);            
                // PRINT("memory::allocate");
                // PRINT(numel);
                int64_t size = allocpos.addr();

                int64_t level = std::max<int64_t>( log2ceil( size ) - log2granularity, 0 );
                        
                // PRINT(level);
                        
                auto node = pool.allocate(level);
                if (!node.first) throw carp::memory::exception(INSUFFICIENT_MEMORY);

                // PRINT(memsize * node.second / sizeof(value_type));

                m_net_allocated += size; // numel;
                m_gross_allocated += (1<<level) * granularity; // / sizeof(MT0);
            
                return granularity * node.second / allocpos.size;
            } // allocate

            void
            release( sizer relpos ) throw ( carp::memory::exception& ) /* override */ {
                // int64_t size = elpos * sizeof(value_type);

                if ( relpos.addr() % granularity != 0 )
                    throw carp::memory::exception(FREEING_UNALLOCATED_MEMORY);
                        
                position_t pos( relpos.addr() / granularity );

                // print_pos( pos, "memory::release::pos"  );
            
                auto released = pool.release(pos);

                m_gross_allocated -= (1<<released) * granularity; // / sizeof(MT0);

                return;
            } // release
        
        }; // class buddy


        // this class manages the continuous segments of the memory
        class local_memory_manager {
            int64_t size;
            int64_t nb_segments;        
            std::vector<int> segments;
        
        public:

            local_memory_manager( int64_t size, int64_t nb_segments, int64_t segment_size )
                : size(size), nb_segments(nb_segments), segments(nb_segments) {
                assert( size >= nb_segments * segment_size );

                for (int q = 0; q < nb_segments; q++ )
                    segments[q] =  (static_cast<float>(size) / nb_segments) * q;
            
            } // cr

            std::vector<int>
            get_segments() {
                return segments;
            } // get_segments

            int operator[] ( int64_t num ) {
                assert(num>0);
                assert(num<nb_segments);
                return segments[num];
            } // operator[]
        
        }; // class volumemanager


        
        class dense : public allocator {
        private:
            int64_t m_pos;
            int64_t m_netallocated;
            int64_t m_grossallocated;
            int64_t m_size;
            
        public:

            dense( sizer allocsize )
                : m_pos(0), m_netallocated(0), m_grossallocated(0), m_size( allocsize.addr() ) { }
                        
            int64_t grossallocated() const /* override */ {
                return m_grossallocated;                
            } // grossallocated
                        
            int64_t netallocated() const /* override */ {
                return m_netallocated;                
            } // netallocated
                        
            int64_t allocate( sizer allocsize ) throw( carp::memory::exception& ) /* override */ {
                int64_t diff = 0;
                if (m_pos % granularity != 0)
                    diff = granularity - m_pos % granularity;
                // PRINT(diff);
                // PRINT(granularity);
                // PRINT(m_pos);                
                // PRINT(granularity - m_pos % granularity);
                
                m_grossallocated += diff;
                m_pos += diff;

                int64_t result = m_pos;                
                                
                int64_t size = allocsize.addr();
//                PRINT(allocsize.addr());                

                if (size + m_pos > m_size) throw carp::memory::exception(INSUFFICIENT_MEMORY);

                m_pos += size;
                m_grossallocated += size;
                m_netallocated += size;
                
                return result / allocsize.size;
            } // allocate
                        
            ~dense() { };
        
        }; // class dense : allocator
        
    } // namespace memory
        
} // namespace CARP






#endif /* __CARP__MEMORY__HPP__ */
// LuM end of file
