angular
  .module('app', ['ui.router','ngMaterial', 'angular-loading-bar','ui.sortable'])
  .controller('IndexController', IndexController)
  .directive('ngType', function () {
    return {
      restrict: 'A',
      link: function (scope, element, attrs) {
        element.attr("type", attrs.ngType);
        console.log('example 1: ' + attrs.ngType);
      }
    }
  }
)
  function IndexController($state){
    this.goto = $state.go
   
  }
  