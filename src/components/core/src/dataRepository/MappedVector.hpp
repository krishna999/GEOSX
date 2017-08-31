/**
 * @file MapVectorContainer.hpp
 * @date Aug 23, 2017
 * @authors settgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_DATAREPOSITORY_MAPVECTORCONTAINER_HPP_
#define SRC_COMPONENTS_CORE_SRC_DATAREPOSITORY_MAPVECTORCONTAINER_HPP_

#include "KeyIndexT.hpp"
#include "SFINAE_Macros.hpp"

/**
 * @class MappedVector
 * This class defines a stl-like container that stores values in an stl vector, and has a map lookup table to
 * access the values by a key. It combines the random access performance of a vector when the index is known,
 * the flexibility of a mapped key lookup O(n) if only the key is known. In addition, a keyIndex can be used for lookup,
 * which will give similar performance to an index lookup after the first use of a keyIndex.
 */
template< typename T, typename T_PTR=T*, typename KEY_TYPE=std::string, typename INDEX_TYPE = int >
class MappedVector
{
public:
  using key_type      = KEY_TYPE ;
  using mapped_type   = T_PTR;

  using LookupMapType          = std::unordered_map<KEY_TYPE, INDEX_TYPE >;
  using value_type             = typename std::pair< KEY_TYPE const, T_PTR >;
  using const_value_type       = typename std::pair< KEY_TYPE const, T const * >;
  using valueContainer         = std::vector<value_type>;
  using const_valueContainer   = std::vector<const_value_type>;
  using pointer                = typename valueContainer::pointer;
  using const_pointer          = typename valueContainer::const_pointer;
  using reference              = typename valueContainer::reference;
  using const_reference        = typename valueContainer::const_reference;
  using size_type              = typename valueContainer::size_type;


  using iterator               = typename valueContainer::iterator;
  using const_iterator         = typename const_valueContainer::const_iterator;
  using reverse_iterator       = typename valueContainer::reverse_iterator;
  using const_reverse_iterator = typename const_valueContainer::const_reverse_iterator;

  using KeyIndex = KeyIndexT<KEY_TYPE const,INDEX_TYPE>;


  MappedVector() = default;
  ~MappedVector() = default;
  MappedVector( MappedVector const & source ) = default ;
  MappedVector & operator=( MappedVector const & source ) = default;
  MappedVector( MappedVector && source ) = default;
  MappedVector & operator=( MappedVector && source ) = default;




  /**
   * @name element access functions
   */
  ///@{


  /**
   * @param index
   * @return pointer to const T
   */
  inline T const * operator[]( INDEX_TYPE index ) const
  {
    return ( index>-1 && index<static_cast<INDEX_TYPE>( m_values.size() ) ) ? const_cast<T const *>(&(*(m_values[index].second))) : nullptr ;
  }

  /**
   *
   * @param index
   * @return pointer to T
   */
  inline T * operator[]( INDEX_TYPE index )
  { return const_cast<T*>( const_cast< MappedVector<T,T_PTR,KEY_TYPE,INDEX_TYPE> const * >(this)->operator[](index) ); }

  /**
   *
   * @param keyName
   * @return pointer to const T
   */
  inline T const * operator[]( KEY_TYPE const & keyName ) const
  {
    typename LookupMapType::const_iterator iter = m_keyLookup.find(keyName);
    return ( iter!=m_keyLookup.end() ? this->operator[]( iter->second ) : nullptr );
  }

  /**
   *
   * @param keyName
   * @return pointer to T
   */
  inline T * operator[]( KEY_TYPE const & keyName )
  { return const_cast<T*>( const_cast< MappedVector<T,T_PTR,KEY_TYPE,INDEX_TYPE> const * >(this)->operator[](keyName) ); }

  /**
   *
   * @param keyIndex
   * @return pointer to const T
   */
  inline T const * operator[]( KeyIndex & keyIndex ) const
  {
    INDEX_TYPE index = keyIndex.Index();

    if( (index==-1) || (m_values[index].first!=keyIndex.Key()) )
    {
      index = getIndex( keyIndex.Key() );
      keyIndex.setIndex(index);
    }

    return this->operator[]( index );
  }

  /**
   *
   * @param keyIndex
   * @return pointer to T
   */
  inline T * operator[]( KeyIndex & keyIndex )
  { return const_cast<T*>( const_cast< MappedVector<T,T_PTR,KEY_TYPE,INDEX_TYPE> const * >(this)->operator[](keyIndex) ); }

  ///@}


  /**
   * @name iterator functions
   */
  ///@{

  /**
   * @return  a read/write iterator that points to the first
   *  element in the in m_objects.
   */
  iterator begin()
  { return m_values.begin(); }

  /**
   * @return  a read-only iterator that points to the first
   *  element in the in m_objects.
   */
  const_iterator begin() const
  { return m_constValues.begin(); }

  /**
   * @return  a read-only iterator that points to the first
   *  element in the in m_objects.
   */
  const_iterator cbegin()
  { return m_constValues.begin(); }

  /**
   * @return  a read/write iterator that points to the last
   *  element in the in m_objects.
   */
  iterator end()
  { return m_values.end(); }

  /**
   * @return  a read-only iterator that points to the last
   *  element in the in m_objects.
   */
  const_iterator end() const
  { return m_constValues.end(); }

  /**
   * @return  a read-only iterator that points to the last
   *  element in the in m_objects.
   */
  const_iterator cend()
  { return m_constValues.end(); }
  ///@}


  /**
   *
   * @param key value of the key to use in the lookup
   * @return index associated with key
   */
  inline INDEX_TYPE getIndex( KEY_TYPE const & key ) const
  {
    typename LookupMapType::const_iterator iter = m_keyLookup.find(key);
    return ( iter!=m_keyLookup.end() ? iter->second : -1 );
  }


  /**
   * @name modifier functions
   */
  ///@{


  T * insert( KEY_TYPE const & keyName, T_PTR source );



  /**
   *  @brief  Remove element at given index
   *  @param  index  index of element to remove.
   *  @return  void
   *
   *  This function will set the element at the given index to nullptr.
   */
  void erase( INDEX_TYPE index )
  {
    m_values[index] = nullptr;
    return;
  }

  /**
   *  @brief  Remove element at given key
   *  @param  key  key of element to remove.
   *  @return  void
   *
   *  This function will set the element at the given key to nullptr.
   */
  void erase( KEY_TYPE const & key )
  {
    typename LookupMapType::const_iterator iter = m_keyLookup.find(key);
    if( iter!=m_keyLookup.end() )
    {
      m_values[iter->second] = nullptr;
    }
    return;
  }

  /**
   *  @brief  Remove element at given key
   *  @param  key  key of element to remove.
   *  @return  void
   *
   *  This function will set the element at the given key to nullptr.
   */
  void erase( KeyIndex & keyIndex )
  {
    INDEX_TYPE index = keyIndex.Index();

    if( (index==-1) || (m_values[index].first!=keyIndex.Key()) )
    {
      index = getIndex( keyIndex.Key() );
      keyIndex.setIndex(index);
    }
    erase( index );
  }

  void clear()
  {
    m_keyLookup.clear();
    m_values.clear();
  }

  ///@}


  inline INDEX_TYPE size() const
  { return m_values.size(); }

  inline const_valueContainer const & values() const
  { return this->m_constValues; }

  inline LookupMapType const & keys() const
  { return m_keyLookup; }


private:
  valueContainer m_values;
  const_valueContainer m_constValues;

  LookupMapType m_keyLookup;
};

template< typename T, typename T_PTR, typename KEY_TYPE, typename INDEX_TYPE >
T * MappedVector<T,T_PTR,KEY_TYPE,INDEX_TYPE>::insert( KEY_TYPE const & keyName , T_PTR source )
{
  typename LookupMapType::iterator iterKeyLookup = m_keyLookup.find(keyName);

  INDEX_TYPE key = -1;
  // if the key was not found, make DataObject<T> and insert
  if( iterKeyLookup == m_keyLookup.end() )
  {
    value_type newEntry = std::make_pair( keyName, std::move( source ) );
    m_values.push_back( std::move( newEntry ) );
    key = m_values.size() - 1;

    m_keyLookup.insert( std::make_pair(keyName,key) );
    m_constValues.push_back( std::make_pair( keyName, &(*(m_values[key].second)) ) );

  }
  // if key was found, make sure it is empty
  else
  {
    key = iterKeyLookup->second;
    if( m_values[key].second==nullptr )
    {
      m_values[key].second = std::move( source );
      m_constValues[key].second =  &(*(m_values[key].second));
    }
    else
    {
      // error?
    }
  }

return &(*(m_values[key].second));
}

#endif /* SRC_COMPONENTS_CORE_SRC_DATAREPOSITORY_MAPVECTORCONTAINER_HPP_ */
