angular
  .module('app')
  .config(routesConfig);

/** @ngInject */
function routesConfig($stateProvider, $urlRouterProvider, $locationProvider) {
  $locationProvider.html5Mode(true).hashPrefix('!');
  $urlRouterProvider.otherwise('/');

  $stateProvider
    .state('app', {
      url: '/',
      component: 'app'
    });
     $stateProvider
    .state('test', {
      url: '/test',
      component: 'test'
    });
    $stateProvider
    .state('reaction', {
      url: '/reaction',
      component: 'reaction'
    });
    $stateProvider
    .state('viewer', {
      url: '/viewer',
      component: 'viewer'
    });

}
