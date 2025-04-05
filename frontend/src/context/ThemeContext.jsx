// src/context/ThemeContext.jsx
import { createContext, useState, useEffect } from 'react';

export const ThemeContext = createContext({
  theme: 'light',
  toggleTheme: () => {},
});

export const ThemeProvider = ({ children }) => {
  // Check if user has a theme preference stored or prefer dark media query
  const getInitialTheme = () => {
    const savedTheme = localStorage.getItem('theme');
    
    if (savedTheme) {
      return savedTheme;
    }
    
    // Check if user prefers dark mode
    return window.matchMedia && 
           window.matchMedia('(prefers-color-scheme: dark)').matches 
           ? 'dark' : 'light';
  };

  const [theme, setTheme] = useState('light'); // Default to light, will be set in useEffect

  // Set the initial theme
  useEffect(() => {
    const initialTheme = getInitialTheme();
    setTheme(initialTheme);
    document.documentElement.setAttribute('data-theme', initialTheme);
  }, []);

  // Theme toggle function
  const toggleTheme = () => {
    const newTheme = theme === 'light' ? 'dark' : 'light';
    setTheme(newTheme);
    document.documentElement.setAttribute('data-theme', newTheme);
    localStorage.setItem('theme', newTheme);
  };

  return (
    <ThemeContext.Provider value={{ theme, toggleTheme }}>
      {children}
    </ThemeContext.Provider>
  );
};