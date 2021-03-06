﻿using System;
using System.Collections.Generic;
using System.Text;
using System.Collections;

namespace Microsoft.Scripting
{
  public class Storage<T> : IAttributesCollection
  {
    #region IAttributesCollection Members

    public void Add(SymbolId name, object value)
    {
      throw new NotImplementedException();
    }

    public bool TryGetValue(SymbolId name, out object value)
    {
      var s = SymbolTable.IdToString(name);
      var fi = Data.GetType().GetField(s);
      if (fi == null)
      {
        value = null;
        return false;
      }
      value = fi.GetValue(Data);
      return true;
    }

    public bool Remove(SymbolId name)
    {
      throw new NotImplementedException();
    }

    public bool ContainsKey(SymbolId name)
    {
      throw new NotImplementedException();
    }

    public object this[SymbolId name]
    {
      get
      {
        throw new NotImplementedException();
      }
      set
      {
        throw new NotImplementedException();
      }
    }

    public IDictionary<SymbolId, object> SymbolAttributes
    {
      get { throw new NotImplementedException(); }
    }

    public void AddObjectKey(object name, object value)
    {
      throw new NotImplementedException();
    }

    public bool TryGetObjectValue(object name, out object value)
    {
      SymbolId n = (SymbolId)name;
      var s = SymbolTable.IdToString(n);
      value = Data.GetType().GetField(s).GetValue(Data);
      return true;
    }

    public bool RemoveObjectKey(object name)
    {
      throw new NotImplementedException();
    }

    public bool ContainsObjectKey(object name)
    {
      throw new NotImplementedException();
    }

    public IDictionary<object, object> AsObjectKeyedDictionary()
    {
      throw new NotImplementedException();
    }

    public int Count
    {
      get { return 0; }
    }

    static readonly ICollection<object> EMPTY = new List<object>();

    public ICollection<object> Keys
    {
      get
      {
        var keys = new List<object>();

        foreach (var fi in Data.GetType().GetFields())
        {
          if (fi.Name != "$parent$")
          {
            keys.Add(fi.Name);
          }
        }

        return keys;
      }
    }

    #endregion

    #region IEnumerable<KeyValuePair<object,object>> Members

    public IEnumerator<KeyValuePair<object, object>> GetEnumerator()
    {
      throw new NotImplementedException();
    }

    #endregion

    #region IEnumerable Members

    IEnumerator IEnumerable.GetEnumerator()
    {
      throw new NotImplementedException();
    }

    #endregion

    public T Data { get; private set; }

    public Storage(T data)
    {
      Data = data;
    }
  }
}
