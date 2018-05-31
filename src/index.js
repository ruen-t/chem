angular
  .module('app', ['ui.router','ngMaterial', 'angular-loading-bar','ui.sortable'])
  .controller('IndexController', IndexController)
  function IndexController($state){
    this.goto = $state.go
   
  }